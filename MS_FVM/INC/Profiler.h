#pragma once

#include <chrono>
#include <iostream>
#include <windows.h>	//compile error is occured when windows.h is not includes before psapi.h
#include <psapi.h>
#include <vector>


class Profiler
{
private:
	static inline PROCESS_MEMORY_COUNTERS_EX memory_recorder_;
	static inline std::vector<size_t> memory_record_;
	static inline std::vector<std::chrono::steady_clock::time_point> time_points_;


public:
	static void set_time_point(void);
	static double get_time_duration(void);

	static void record_Consumed_Memory(void);
	static void print_Consumed_Memory(void);
	static void record_Consumed_Memory_And_Time(void);
	static void print_Consumed_Memory_And_Time(void);	
};


#define SET_TIME_POINT Profiler::set_time_point()
#define GET_TIME_DURATION Profiler::get_time_duration()
#define RECORD_CONSUMED_MEMORY Profiler::record_Consumed_Memory()
#define PRINT_CONSUMED_MEMORY Profiler::print_Consumed_Memory()
#define RECORD_CONSUMED_MEMORY_AND_TIME Profiler::record_Consumed_Memory_And_Time()
#define PRINT_CONSUMED_MEMORY_AND_TIME Profiler::print_Consumed_Memory_And_Time()
