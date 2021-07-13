
#include "../INC/Profiler.h"

void Profiler::record_Consumed_Memory(void){		
	GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&memory_recorder_, sizeof(memory_recorder_));

	memory_record_.push_back(memory_recorder_.PrivateUsage);
}

void Profiler::print_Consumed_Memory(void){
	GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&memory_recorder_, sizeof(memory_recorder_));

	const auto consumed_memory = static_cast<size_t>((memory_recorder_.PrivateUsage - memory_record_.back()) / 1024.0 / 1024.0);
	std::cout << "\tconsummed_memory : " << consumed_memory << " MB\n";

	memory_record_.pop_back();
}

void Profiler::set_time_point(void){
	time_points_.push_back(std::chrono::steady_clock::now());
}

double Profiler::get_time_duration(void){
	std::chrono::duration<double> time_duration = std::chrono::steady_clock::now() - time_points_.back();	
	time_points_.pop_back();
	return time_duration.count();
}

void Profiler::record_Consumed_Memory_And_Time(void){
	record_Consumed_Memory();
	set_time_point();
}

void Profiler::print_Consumed_Memory_And_Time(void) {
	if (memory_record_.size() == 0 || time_points_.size() == 0)
		throw std::runtime_error("there are no record to print!");

	GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&memory_recorder_, sizeof(memory_recorder_));

	const auto consummed_memory = static_cast<size_t>((memory_recorder_.PrivateUsage - memory_record_.back()) / 1024.0 / 1024.0);
	const std::chrono::duration<double> consummed_time = std::chrono::steady_clock::now() - time_points_.back();

	std::cout << "\tconsummed memory : " << consummed_memory << " MB" << "  \tconsummed time : " << consummed_time.count() << " s\n";

	memory_record_.pop_back();
	time_points_.pop_back();
}