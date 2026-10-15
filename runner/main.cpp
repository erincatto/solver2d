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
	"jacobi", "pgs", "pgs_ngs", "pgs_ngs_block", "pgs_soft", "soft_step", "tgs_sticky", "tgs_soft", "tgs_ngs", "xpbd", "dvi",
};

static_assert(S2_ARRAY_COUNT(s_solverNames) == s2_solverTypeCount);

struct Options
{
	std::vector<int> samples;
	bool solvers[s2_solverTypeCount] = {};
	int stepCount = 600;
	float hertz = 60.0f;
	int primaryIterations = 4;
	int secondaryIterations = 2;
	bool warmStart = true;
	bool dumpBodies = false;
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
	printf("  -iters=<n>            primary iterations, substeps for substepping solvers (default 4)\n");
	printf("  -extra=<n>            secondary iterations, relax or position (default 2)\n");
	printf("  -nowarm               disable warm starting\n");
	printf("  -bodies               also dump per body state each step\n");
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

	if (options.samples.empty() || options.stepCount <= 0 || options.hertz <= 0.0f)
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
					  "max_joint_error,mean_joint_error,max_limit_error,kinetic,potential,energy,max_speed,max_angular_speed\n");

	Sample* sample = g_sampleEntries[sampleIndex].createFcn(settings, solverType);
	Probe* probe = CreateProbe();
	std::vector<ProbeBody> bodies;

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

		fprintf(stepFile, "%d,%.6g,%.4f,%d,%d,%d,%.6g,%.6g,%.6g,%d,%.6g,%.6g,%.6g,%.8g,%.8g,%.8g,%.6g,%.6g\n", step, time, ms,
				s.bodyCount, s.contactCount, s.pointCount, s.maxPenetration, s.meanPenetration, s.normalImpulse, s.jointCount,
				s.maxJointError, s.meanJointError, s.maxLimitError, s.kineticEnergy, s.potentialEnergy, energy, s.maxSpeed,
				s.maxAngularSpeed);

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
						 "final_joint_error,max_limit_error,initial_energy,final_energy,final_kinetic,total_ms\n");

	printf("steps=%d hertz=%g iters=%d extra=%d warm=%d out=%s\n", options.stepCount, options.hertz, options.primaryIterations,
		   options.secondaryIterations, options.warmStart ? 1 : 0, outDir.string().c_str());

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

			fprintf(summaryFile, "\"%s\",%s,%d,%g,%d,%d,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.8g,%.8g,%.8g,%.3f\n", full.c_str(),
					s_solverNames[i], options.stepCount, options.hertz, options.primaryIterations, options.secondaryIterations,
					options.warmStart ? 1 : 0, r.maxPenetration, r.finalPenetration, r.maxJointError, r.finalJointError,
					r.maxLimitError, r.initialEnergy, r.finalEnergy, r.finalKinetic, r.totalMs);
		}
	}

	fclose(summaryFile);
	return failures == 0 ? 0 : 1;
}
