#pragma once
#include "Semi_Discrete_Equation.h"
#include "Time_Integral_Method.h"
#include "Solve_Condition.h"

template <typename T>
class Discrete_Equation
{
    static_require(ms::is_Time_Intg_Method<T>, "Wrong time integral method");

    using Time_Integ_Method = T;

public:
 
    template<typename End_Cond, typename Post_Cond, typename SDE>
    static void solve(const SDE& semi_discrete_eq, std::vector<typename SDE::Solution>& solutions) {
        static_require(ms::is_End_Condtion<End_Cond>, "Wrong end condition");
        static_require(ms::is_Post_Condtion<Post_Cond>, "Wrong post condition");

        // post initial solution
        while (true) {
            auto time_step = semi_discrete_eq.calculate_time_step(solutions);

            if (End_Cond::check(time_step)) {
                Time_Integ_Method::update_solutions(semi_discrete_eq, solutions, time_step);
                //post last solution
                break;
            }

            if (Post_Cond::check(time_step)) {
                Time_Integ_Method::update_solutions(semi_discrete_eq, solutions, time_step);
                //post
            }
            else
                Time_Integ_Method::update_solutions(semi_discrete_eq, solutions, time_step);
        }
    }
};