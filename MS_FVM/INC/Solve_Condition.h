#pragma once
#include "Log.h"

#include <type_traits>


class SEC {}; // Solve End Condition

template<double target_time>
class End_By_Time : public SEC
{
public:
    static bool check(const double current_time, const double time_step) {
        const double expect_time = current_time + time_step;        
        if (target_time <= expect_time) {
            Log::content_ << "current time: " << std::to_string(target_time) + "s  (100.00%)";
            return true;
        }
        else {
            Log::content_ << "current time: " << std::to_string(expect_time) + "s  ";
            Log::content_ << std::fixed << std::setprecision(3) << "(" << expect_time * 100 / target_time << "%) " << std::defaultfloat << std::setprecision(6);
            return false;
        }
    }

    static void adjust(double& current_time, double& time_step) {
        const auto expect_time = current_time + time_step;
        const auto exceed_time = expect_time - target_time;

        time_step -= exceed_time;
        current_time = target_time;        
    }
};


class SPC {};   // Solve Post Condition

template<double post_time_step>
class Post_By_Time : public SPC
{
public:
    static bool check(const double current_time, const double time_step) {
        const auto target_time = num_post_ * post_time_step;
        const auto expect_time = current_time + time_step;
        if (target_time <= expect_time)
            return true;
        else
            return false;
    }

    static void adjust(double& current_time, double& time_step) {
        const auto expect_time = current_time + time_step;
        const auto target_time = num_post_ * post_time_step;
        const auto exceed_time = expect_time - target_time;

        time_step -= exceed_time;
        current_time = target_time;

        num_post_++;
    }

private:
    inline static size_t num_post_ = 1;
};

namespace ms {
    template <typename T>
    inline constexpr bool is_solve_end_condtion = std::is_base_of_v<SEC, T>;
    template <typename T>
    inline constexpr bool is_solve_post_condtion = std::is_base_of_v<SPC, T>;

}
