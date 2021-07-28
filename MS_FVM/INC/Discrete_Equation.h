#pragma once
#include "Semi_Discrete_Equation.h"
#include "Time_Integral_Method.h"
#include "Time_Step_Method.h"
#include "Solve_Condition.h"
#include "Post.h"

template <typename Time_Integral_Method>
class Discrete_Equation
{
    static_require(ms::is_time_integral_method<Time_Integral_Method>, "It should be time integral method");

public: 
    template<typename Time_Step_Method, typename Solve_End_Condition, typename Solve_Post_Condition, typename Post, typename SDE, typename Solution>
    static void solve(const SDE& semi_discrete_eq, std::vector<Solution>& solutions) {
        static_require(ms::is_solve_end_condtion<Solve_End_Condition>,      "It should be solve end condition");
        static_require(ms::is_solve_post_condtion<Solve_Post_Condition>,    "It should be solve post condition");
         
        double current_time = 0.0;
        Post::syncronize_time(current_time);
        Post::solution(solutions, "initial");

        Log::content_ << "================================================================================\n";
        Log::content_ << "\t\t\t\t Solving\n";
        Log::content_ << "================================================================================\n\t\t\t\t\t\t";

        SET_TIME_POINT;
        while (true) {
            SET_TIME_POINT;
            auto time_step = semi_discrete_eq.calculate_time_step<Time_Step_Method>(solutions);
             
            if (Solve_End_Condition::inspect(current_time, time_step)) {
                Time_Integral_Method::update_solutions(semi_discrete_eq, solutions, time_step);
                current_time += time_step;
                Log::content_ << "time/update: " << std::to_string(GET_TIME_DURATION) << "s   \t";

                Log::content_ << "current time: " << std::to_string(current_time) + "s  (100.00%)\n";
                Post::solution(solutions, "final");
                break;
            }

            if (Solve_Post_Condition::inspect(current_time, time_step)) {
                Time_Integral_Method::update_solutions(semi_discrete_eq, solutions, time_step);
                current_time += time_step;

                Post::solution(solutions);
            }
            else {
                Time_Integral_Method::update_solutions(semi_discrete_eq, solutions, time_step);
                current_time += time_step;
            }          

            Log::content_ << "time/update: " << std::to_string(GET_TIME_DURATION) << "s   \t";
            Log::print();                
        }

        Log::content_ << "================================================================================\n";
        Log::content_ << "\t\t\t Total ellapsed time: " << GET_TIME_DURATION << "s\n";
        Log::content_ << "================================================================================\n\n";
        Log::print();
    }
};