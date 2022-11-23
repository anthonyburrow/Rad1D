#include <iostream>
#include <catch2/catch_test_macros.hpp>

#include "test_main.hpp"
#include "../../src/blackbody.hpp"

using namespace std;

namespace myLib
{
    void testAll()
    {
        TEST_CASE( "Planck function tested", "[planck]" ) {
            REQUIRE( planck(5000.0, 6000.0) == 1.0 );
        }
    }
}

