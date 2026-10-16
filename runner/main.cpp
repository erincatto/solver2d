// SPDX-FileCopyrightText: 2024 Erin Catto
// SPDX-License-Identifier: MIT

// Headless scene runner. Steps sample scenes for each solver and writes per step CSV.

#include "probe.h"
#include "sample.h"
#include "settings.h"

#include "solver2d/math.h"

#include <chrono>
#include <ctype.h>
#include <filesystem>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

static const char* s_solverNames[s2_solverTypeCount] = {
	"jacobi", "pgs", "pgs_ngs", "pgs_ngs_block", "pgs_soft", "soft_step", "tgs_sticky", "tgs_soft", "tgs_ngs", "xpbd", "dvi", "soft_step_d", "soft_step_seed",
};

static_assert(S2_ARRAY_COUNT(s_solverNames) == s2_solverTypeCount);

struct Options
{
	std::vector<int> samples;
	bool solvers[s2_solverTypeCount] = {};
	int stepCount = 600;
	float hertz = 60.0f;
	float jointHertz = 60.0f;
	float jointHertzD = 480.0f;
	float jointHertzCapD = 2.0f;
	int jointCRIterationsD = 0;
	bool directRelaxD = true;
	float jointHertzSeed = 60.0f;
	float jointHertzCapSeed = 0.25f;
	int seedMode = 0;
	int primaryIterations = 4;
	int secondaryIterations = 2;
	bool warmStart = true;
	bool dumpBodies = false;
	bool measureCR = false;
	std::string outDir = "build/runner_out";
};

static std::string Lower(const char* s)
{
	std::string r = s;
	for (char& c : r)
	{
		c = (char)tolower((unsigned char)c);
	}
	return r;
}

static std::string FullName(int index)
{
	return std::string(g_sampleEntries[index].category) + "/" + g_sampleEntries[index].name;
}

static std::string FileName(const std::string& s)
{
	std::string r;
	for (char c : s)
	{
		if (isalnum((unsigned char)c))
		{
			r += (char)tolower((unsigned char)c);
		}
		else if (r.empty() == false && r.back() != '_')
		{
			r += '_';
		}
	}

	while (r.empty() == false && r.back() == '_')
	{
		r.pop_back();
	}

	return r;
}

static void PrintUsage()
{
	printf("usage: runner -sample=<index|name|all> [options]\n");
	printf("  -list                 list samples and solvers\n");
	printf("  -sample=<s>           sample index, full or partial name (case insensitive), or all\n");
	printf("  -solvers=<a,b,...>    solver names or all (default all)\n");
	printf("  -steps=<n>            steps to run (default 600)\n");
	printf("  -hertz=<f>            step rate (default 60)\n");
	printf("  -joint_hertz=<f>      Soft Step joint stiffness, capped at a quarter of the substep rate (default 60)\n");
	printf("  -joint_hertz_d=<f>    Soft Step D joint stiffness (default 480)\n");
	printf("  -joint_cap_d=<f>      Soft Step D joint hertz cap as a fraction of the substep rate, 0 for none (default 2)\n");
	printf("  -cr_iters_d=<n>       Soft Step D solves joints by n conjugate residual iterations, 0 for exact (default 0)\n");
	printf("  -gs_relax_d           Soft Step D relaxes joints with Soft Step's Gauss-Seidel pass, not the direct solve\n");
	printf("  -joint_hertz_seed=<f> Soft Step Seed joint stiffness (default 60)\n");
	printf("  -joint_cap_seed=<f>   Soft Step Seed joint hertz cap as a fraction of the substep rate, 0 for none (default 0.25)\n");
	printf("  -seed=<s>             Soft Step Seed exact solve: first (substep), every (substep), or step (default first)\n");
	printf("  -iters=<n>            primary iterations, substeps for substepping solvers (default 4)\n");
	printf("  -extra=<n>            secondary iterations, relax or position (default 2)\n");
	printf("  -nowarm               disable warm starting\n");
	printf("  -bodies               also dump per body state each step\n");
	printf("  -cr                   Soft Step D measures conjugate residual against its exact joint solves (slow)\n");
	printf("  -out=<dir>            output directory (default build/runner_out)\n");
}

static void PrintList()
{
	printf("samples:\n");
	for (int i = 0; i < g_sampleCount; ++i)
	{
		printf("  %2d  %s\n", i, FullName(i).c_str());
	}

	printf("solvers:\n");
	for (int i = 0; i < s2_solverTypeCount; ++i)
	{
		printf("  %s\n", s_solverNames[i]);
	}
}

static bool ParseSample(const char* arg, std::vector<int>& samples)
{
	std::string key = Lower(arg);
	if (key == "all")
	{
		for (int i = 0; i < g_sampleCount; ++i)
		{
			samples.push_back(i);
		}
		return true;
	}

	char* end;
	long index = strtol(arg, &end, 10);
	if (*end == 0 && end != arg)
	{
		if (index < 0 || index >= g_sampleCount)
		{
			fprintf(stderr, "sample index %ld out of range\n", index);
			return false;
		}
		samples.push_back(int(index));
		return true;
	}

	// Exact match wins so "Joints/Ragdoll" is not ambiguous with "Ragdoll Stress"
	std::vector<int> partial;
	for (int i = 0; i < g_sampleCount; ++i)
	{
		std::string full = Lower(FullName(i).c_str());
		std::string name = Lower(g_sampleEntries[i].name);
		if (key == full || key == name)
		{
			samples.push_back(i);
			return true;
		}

		if (full.find(key) != std::string::npos)
		{
			partial.push_back(i);
		}
	}

	if (partial.size() == 1)
	{
		samples.push_back(partial[0]);
		return true;
	}

	if (partial.empty())
	{
		fprintf(stderr, "no sample matches '%s'\n", arg);
	}
	else
	{
		fprintf(stderr, "'%s' is ambiguous:\n", arg);
		for (int i : partial)
		{
			fprintf(stderr, "  %2d  %s\n", i, FullName(i).c_str());
		}
	}

	return false;
}

static bool ParseSolvers(const char* arg, bool* solvers)
{
	std::string list = Lower(arg);
	size_t start = 0;
	while (start <= list.size())
	{
		size_t comma = list.find(',', start);
		std::string token = list.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
		start = comma == std::string::npos ? list.size() + 1 : comma + 1;

		if (token.empty())
		{
			continue;
		}

		if (token == "all")
		{
			for (int i = 0; i < s2_solverTypeCount; ++i)
			{
				solvers[i] = true;
			}
			continue;
		}

		int found = -1;
		for (int i = 0; i < s2_solverTypeCount; ++i)
		{
			if (token == s_solverNames[i])
			{
				found = i;
				break;
			}
		}

		if (found < 0)
		{
			fprintf(stderr, "unknown solver '%s'\n", token.c_str());
			return false;
		}

		solvers[found] = true;
	}

	return true;
}

static bool ParseArgs(int argc, char** argv, Options& options)
{
	bool anySolver = false;
	for (int i = 1; i < argc; ++i)
	{
		const char* arg = argv[i];
		const char* value = strchr(arg, '=');
		value = value != nullptr ? value + 1 : "";

		if (strcmp(arg, "-list") == 0)
		{
			PrintList();
			exit(0);
		}
		else if (strncmp(arg, "-sample=", 8) == 0)
		{
			if (ParseSample(value, options.samples) == false)
			{
				return false;
			}
		}
		else if (strncmp(arg, "-solvers=", 9) == 0)
		{
			if (ParseSolvers(value, options.solvers) == false)
			{
				return false;
			}
			anySolver = true;
		}
		else if (strncmp(arg, "-steps=", 7) == 0)
		{
			options.stepCount = atoi(value);
		}
		else if (strncmp(arg, "-hertz=", 7) == 0)
		{
			options.hertz = (float)atof(value);
		}
		else if (strncmp(arg, "-joint_hertz=", 13) == 0)
		{
			options.jointHertz = (float)atof(value);
		}
		else if (strncmp(arg, "-joint_hertz_d=", 15) == 0)
		{
			options.jointHertzD = (float)atof(value);
		}
		else if (strncmp(arg, "-joint_cap_d=", 13) == 0)
		{
			options.jointHertzCapD = (float)atof(value);
		}
		else if (strncmp(arg, "-cr_iters_d=", 12) == 0)
		{
			options.jointCRIterationsD = atoi(value);
		}
		else if (strcmp(arg, "-gs_relax_d") == 0)
		{
			options.directRelaxD = false;
		}
		else if (strncmp(arg, "-joint_hertz_seed=", 18) == 0)
		{
			options.jointHertzSeed = (float)atof(value);
		}
		else if (strncmp(arg, "-joint_cap_seed=", 16) == 0)
		{
			options.jointHertzCapSeed = (float)atof(value);
		}
		else if (strncmp(arg, "-seed=", 6) == 0)
		{
			const char* modes[] = {"first", "every", "step"};
			static_assert(S2_ARRAY_COUNT(modes) == s2_seedModeCount);
			options.seedMode = -1;
			for (int i = 0; i < s2_seedModeCount; ++i)
			{
				if (strcmp(value, modes[i]) == 0)
				{
					options.seedMode = i;
				}
			}

			if (options.seedMode < 0)
			{
				fprintf(stderr, "unknown seed '%s'\n", value);
				return false;
			}
		}
		else if (strncmp(arg, "-iters=", 7) == 0)
		{
			options.primaryIterations = atoi(value);
		}
		else if (strncmp(arg, "-extra=", 7) == 0)
		{
			options.secondaryIterations = atoi(value);
		}
		else if (strcmp(arg, "-nowarm") == 0)
		{
			options.warmStart = false;
		}
		else if (strcmp(arg, "-bodies") == 0)
		{
			options.dumpBodies = true;
		}
		else if (strcmp(arg, "-cr") == 0)
		{
			options.measureCR = true;
		}
		else if (strncmp(arg, "-out=", 5) == 0)
		{
			options.outDir = value;
		}
		else
		{
			fprintf(stderr, "unknown argument '%s'\n", arg);
			return false;
		}
	}

	if (anySolver == false)
	{
		ParseSolvers("all", options.solvers);
	}

	if (options.samples.empty() || options.stepCount <= 0 || options.hertz <= 0.0f || options.jointHertz <= 0.0f ||
		options.jointHertzD <= 0.0f || options.jointHertzCapD < 0.0f || options.jointCRIterationsD < 0 ||
		options.jointHertzSeed <= 0.0f || options.jointHertzCapSeed < 0.0f)
	{
		return false;
	}

	return true;
}

struct RunSummary
{
	float maxPenetration = 0.0f;
	float finalPenetration = 0.0f;
	float maxJointError = 0.0f;
	float finalJointError = 0.0f;
	float maxLimitError = 0.0f;
	float initialEnergy = 0.0f;
	float finalEnergy = 0.0f;
	float finalKinetic = 0.0f;
	double totalMs = 0.0;
};

static bool RunSolver(const Options& options, int sampleIndex, s2SolverType solverType, const std::filesystem::path& dir,
					  RunSummary& summary)
{
	Settings settings;
	settings.hertz = options.hertz;
	settings.timeStep = 1.0f / options.hertz;
	settings.jointHertz = options.jointHertz;
	settings.jointHertzD = options.jointHertzD;
	settings.jointHertzCapD = options.jointHertzCapD;
	settings.jointCRIterationsD = options.jointCRIterationsD;
	settings.directRelaxD = options.directRelaxD;
	settings.jointHertzSeed = options.jointHertzSeed;
	settings.jointHertzCapSeed = options.jointHertzCapSeed;
	settings.seedMode = options.seedMode;
	settings.primaryIterations = options.primaryIterations;
	settings.secondaryIterations = options.secondaryIterations;
	settings.enableWarmStarting = options.warmStart;

	// Keeps Sample::Step from issuing debug draw
	for (bool& enabled : settings.enabledSolvers)
	{
		enabled = false;
	}

	std::string stepPath = (dir / (std::string(s_solverNames[solverType]) + ".csv")).string();
	FILE* stepFile = fopen(stepPath.c_str(), "w");
	if (stepFile == nullptr)
	{
		fprintf(stderr, "cannot write %s\n", stepPath.c_str());
		return false;
	}

	FILE* bodyFile = nullptr;
	if (options.dumpBodies)
	{
		std::string bodyPath = (dir / (std::string(s_solverNames[solverType]) + "_bodies.csv")).string();
		bodyFile = fopen(bodyPath.c_str(), "w");
		if (bodyFile == nullptr)
		{
			fprintf(stderr, "cannot write %s\n", bodyPath.c_str());
			fclose(stepFile);
			return false;
		}
		fprintf(bodyFile, "step,time,body,x,y,angle,vx,vy,w\n");
	}

	fprintf(stepFile, "step,time,ms,bodies,contacts,points,max_penetration,mean_penetration,normal_impulse,joints,"
					  "max_joint_error,mean_joint_error,max_limit_error,kinetic,potential,energy,max_speed,max_angular_speed,"
					  "max_angle_error%s\n",
			options.measureCR ? ",cr_soft_solves,cr_soft_it1,cr_soft_it2,cr_soft_e8,cr_soft_e32,cr_soft_exact,cr_relax_solves,"
								"cr_relax_it1,cr_relax_it2,cr_relax_e8,cr_relax_e32,cr_relax_exact"
							  : "");

	Sample* sample = g_sampleEntries[sampleIndex].createFcn(settings, solverType);
	Probe* probe = CreateProbe();
	std::vector<ProbeBody> bodies;

	if (options.measureCR)
	{
		Probe_EnableJointCR(sample->m_worldId);
	}

	// Step 0 is the initial state
	for (int step = 0; step <= options.stepCount; ++step)
	{
		Probe_Snapshot(probe, sample->m_worldId);

		double ms = 0.0;
		if (step > 0)
		{
			auto t1 = std::chrono::high_resolution_clock::now();
			sample->Step(settings, s2Color{1.0f, 1.0f, 1.0f, 1.0f});
			auto t2 = std::chrono::high_resolution_clock::now();
			ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
		}

		ProbeStep s = Probe_Measure(probe, sample->m_worldId);
		float time = step * settings.timeStep;
		float energy = s.kineticEnergy + s.potentialEnergy;

		fprintf(stepFile, "%d,%.6g,%.4f,%d,%d,%d,%.6g,%.6g,%.6g,%d,%.6g,%.6g,%.6g,%.8g,%.8g,%.8g,%.6g,%.6g,%.6g", step, time, ms,
				s.bodyCount, s.contactCount, s.pointCount, s.maxPenetration, s.meanPenetration, s.normalImpulse, s.jointCount,
				s.maxJointError, s.meanJointError, s.maxLimitError, s.kineticEnergy, s.potentialEnergy, energy, s.maxSpeed,
				s.maxAngularSpeed, s.maxAngleError);

		if (options.measureCR)
		{
			fprintf(stepFile, ",%d,%d,%d,%.4g,%.4g,%.4g,%d,%d,%d,%.4g,%.4g,%.4g", s.crSoftSolves, s.crSoftIterations1,
					s.crSoftIterations2, s.crSoftError8, s.crSoftError32, s.crSoftExactResidual, s.crRelaxSolves,
					s.crRelaxIterations1, s.crRelaxIterations2, s.crRelaxError8, s.crRelaxError32, s.crRelaxExactResidual);
		}

		fprintf(stepFile, "\n");

		if (bodyFile != nullptr)
		{
			int count = Probe_GetBodies(sample->m_worldId, bodies.data(), int(bodies.size()));
			if (count > int(bodies.size()))
			{
				bodies.resize(count);
				count = Probe_GetBodies(sample->m_worldId, bodies.data(), int(bodies.size()));
			}

			for (int i = 0; i < count; ++i)
			{
				const ProbeBody& b = bodies[i];
				fprintf(bodyFile, "%d,%.6g,%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n", step, time, b.index, b.position.x, b.position.y,
						b.angle, b.linearVelocity.x, b.linearVelocity.y, b.angularVelocity);
			}
		}

		if (step == 0)
		{
			summary.initialEnergy = energy;
		}

		summary.maxPenetration = s.maxPenetration > summary.maxPenetration ? s.maxPenetration : summary.maxPenetration;
		summary.maxJointError = s.maxJointError > summary.maxJointError ? s.maxJointError : summary.maxJointError;
		summary.maxLimitError = s.maxLimitError > summary.maxLimitError ? s.maxLimitError : summary.maxLimitError;
		summary.finalPenetration = s.maxPenetration;
		summary.finalJointError = s.maxJointError;
		summary.finalEnergy = energy;
		summary.finalKinetic = s.kineticEnergy;
		summary.totalMs += ms;
	}

	DestroyProbe(probe);
	delete sample;

	fclose(stepFile);
	if (bodyFile != nullptr)
	{
		fclose(bodyFile);
	}

	return true;
}

int main(int argc, char** argv)
{
	Options options;
	if (ParseArgs(argc, argv, options) == false)
	{
		PrintUsage();
		return 1;
	}

	std::filesystem::path outDir = options.outDir;
	std::error_code ec;
	std::filesystem::create_directories(outDir, ec);

	std::string summaryPath = (outDir / "summary.csv").string();
	FILE* summaryFile = fopen(summaryPath.c_str(), "w");
	if (summaryFile == nullptr)
	{
		fprintf(stderr, "cannot write %s\n", summaryPath.c_str());
		return 1;
	}

	fprintf(summaryFile, "sample,solver,steps,hertz,iters,extra,warm_start,max_penetration,final_penetration,max_joint_error,"
						 "final_joint_error,max_limit_error,initial_energy,final_energy,final_kinetic,total_ms,joint_hertz,joint_hertz_d,"
						 "joint_cap_d,cr_iters_d,joint_hertz_seed,joint_cap_seed,seed\n");

	printf("steps=%d hertz=%g iters=%d extra=%d warm=%d joint_hertz=%g joint_hertz_d=%g joint_cap_d=%g cr_iters_d=%d out=%s\n",
		   options.stepCount, options.hertz, options.primaryIterations, options.secondaryIterations, options.warmStart ? 1 : 0,
		   options.jointHertz, options.jointHertzD, options.jointHertzCapD, options.jointCRIterationsD, outDir.string().c_str());

	int failures = 0;
	for (int sampleIndex : options.samples)
	{
		std::string full = FullName(sampleIndex);
		std::filesystem::path dir = outDir / FileName(full);
		std::filesystem::create_directories(dir, ec);

		printf("\n%s\n", full.c_str());
		printf("  %-14s %10s %10s %10s %10s %12s %12s %9s\n", "solver", "max_pen", "final_pen", "max_jerr", "final_jerr",
			   "energy0", "energy", "ms");

		for (int i = 0; i < s2_solverTypeCount; ++i)
		{
			if (options.solvers[i] == false)
			{
				continue;
			}

			RunSummary r;
			if (RunSolver(options, sampleIndex, s2SolverType(i), dir, r) == false)
			{
				failures += 1;
				continue;
			}

			printf("  %-14s %10.4g %10.4g %10.4g %10.4g %12.6g %12.6g %9.1f\n", s_solverNames[i], r.maxPenetration,
				   r.finalPenetration, r.maxJointError, r.finalJointError, r.initialEnergy, r.finalEnergy, r.totalMs);

			fprintf(summaryFile, "\"%s\",%s,%d,%g,%d,%d,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.8g,%.8g,%.8g,%.3f,%g,%g,%g,%d,%d,%g,%g,%d\n", full.c_str(),
					s_solverNames[i], options.stepCount, options.hertz, options.primaryIterations, options.secondaryIterations,
					options.warmStart ? 1 : 0, r.maxPenetration, r.finalPenetration, r.maxJointError, r.finalJointError,
					r.maxLimitError, r.initialEnergy, r.finalEnergy, r.finalKinetic, r.totalMs, options.jointHertz,
					options.jointHertzD, options.jointHertzCapD, options.jointCRIterationsD, options.directRelaxD ? 1 : 0,
					options.jointHertzSeed, options.jointHertzCapSeed, options.seedMode);
		}
	}

	fclose(summaryFile);
	return failures == 0 ? 0 : 1;
}
