#pragma once
#include "Semi_Discrete_Equation.h"
#include "Time_Integral_Method.h"
#include "Time_Step_Method.h"
#include "Solve_Condition.h"
#include "Post.h"

template <typename TIM>
class Discrete_Equation
{
    static_require(ms::is_time_integral_method<TIM>, "Wrong time integral method");

    using Time_Integral_Method = TIM;

public: 
    template<typename Time_Step_Method, typename Solve_End_Condition, typename Solve_Post_Condition, typename SDE, typename Solution>
    static void solve(const SDE& semi_discrete_eq, std::vector<Solution>& solutions) {
        static_require(ms::is_solve_end_condtion<Solve_End_Condition>, "Wrong solve end condition");
        static_require(ms::is_solve_post_condtion<Solve_Post_Condition>, "Wrong solve post condition");
         
        double current_time = 0.0;
        Post::solution(solutions, current_time, "initial");

        Log::content_ << "================================================================================\n";
        Log::content_ << "\t\t\t\t Solving\n";
        Log::content_ << "================================================================================\n";

        SET_TIME_POINT;
        while (true) {
            SET_TIME_POINT;
            auto time_step = semi_discrete_eq.calculate_time_step<Time_Step_Method>(solutions);
             
            if (Solve_End_Condition::check(current_time, time_step)) {
                Solve_End_Condition::adjust(current_time, time_step);
                Time_Integral_Method::update_solutions(semi_discrete_eq, solutions, time_step);
                Log::content_ << std::left << "\t time/iter: " << GET_TIME_DURATION << "s\n";
                Post::solution(solutions, current_time, "final");
                break;
            }

            if (Solve_Post_Condition::check(current_time, time_step)) {
                Solve_Post_Condition::adjust(current_time, time_step);
                Time_Integral_Method::update_solutions(semi_discrete_eq, solutions, time_step);
                Post::solution(solutions, current_time);
            }
            else {
                Time_Integral_Method::update_solutions(semi_discrete_eq, solutions, time_step);
                current_time += time_step;
            }          

            Log::content_ << std::left << "\t time/iter: " << GET_TIME_DURATION << "s\n";
            Log::print();
                
        }

        Log::content_ << "================================================================================\n";
        Log::content_ << "\t\t\t Total ellapsed time: " << GET_TIME_DURATION << "s\n";
        Log::content_ << "================================================================================\n\n";
        Log::print();
    }
};