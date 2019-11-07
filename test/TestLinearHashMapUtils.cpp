//
// Created by xetql on 11/6/19.
//

#include <gtest/gtest.h>
#include "Utils.hpp"

namespace
{

TEST(LinearHashMap, isSearchable) {
    LinearHashMap <int, int, 8> hashMap = { std::make_pair(0, 101)
            ,std::make_pair(1, 202)
            ,std::make_pair(2, 303)
            ,std::make_pair(3, 404)
            ,std::make_pair(4, 505)
            ,std::make_pair(5, 606)
            ,std::make_pair(6, 707)
            ,std::make_pair(7, 808)};
    EXPECT_EQ(303, (*search_in_linear_hashmap<int, int, 8>(hashMap, 2)).second);
    EXPECT_EQ(hashMap.end(), (search_in_linear_hashmap<int, int, 8>(hashMap, 10)));
}

TEST(LinearHashMapMutableRef, isMutable) {
    LinearHashMap <int, bool, 8> hashMap = { std::make_pair(0, true)
                                            ,std::make_pair(1, true)
                                            ,std::make_pair(2, true)
                                            ,std::make_pair(3, true)
                                            ,std::make_pair(4, true)
                                            ,std::make_pair(5, true)
                                            ,std::make_pair(6, true)
                                            ,std::make_pair(7, true)};
    bool& val = (*search_in_linear_hashmap<int, bool, 8>(hashMap, 2)).second;
    val = false;
    EXPECT_EQ(false, (*search_in_linear_hashmap<int, bool, 8>(hashMap, 2)).second);
}

} //end of namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}