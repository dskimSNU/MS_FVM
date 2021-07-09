#pragma once
#include <type_traits>

class End_Condition {};

template<double target_time>
class End_By_Time : public End_Condition
{
public:
    static bool check(double& time_step) {
        static double time = 0.0;
        time += time_step;
        if (target_time <= time) {
            const auto exceed_time = time - target_time;
            time_step -= exceed_time;
            return true;
        }
        return false;
    }
};


class Post_Condition {};

template<double target_time>
class Post_By_Time : public Post_Condition
{
public:
    static bool check(double& time_step) {
        static double time = 0.0;
        time += time_step;
        if (target_time <= time) {
            const auto exceed_time = time - target_time;
            time_step -= exceed_time;
            return true;
        }
        return false;
    }
};

namespace ms {
    template <typename T>
    inline constexpr bool is_End_Condtion = std::is_base_of_v<End_Condition, T>;
    template <typename T>
    inline constexpr bool is_Post_Condtion = std::is_base_of_v<Post_Condition, T>;

}
