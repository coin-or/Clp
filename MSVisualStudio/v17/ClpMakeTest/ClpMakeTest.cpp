// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file ClpMakeTest.cpp
    \brief Main routine for running quick tests that need no tools beyond the build environment, similar to 'make test'.
*/

#pragma once

#include "CoinError.hpp" // for WindowsErrorPopupBlocker()
#include <cstdio> // for printf

// The order below of including the test cpp files is important to prevent 
// ClpSimplex from being defined multiple times.

#define main main_test_racing_lp
#include "..\test\test_racing_lp.cpp"
#undef main

#define main main_test_barrier_control
#include "..\test\test_barrier_control.cpp"
#undef main

// Run the 'test' suite in order, aborting on the first failure.
int test(int argc, char** argv)
{
    int returnCode = 0;

    //"Exit code: 0 = all checks passed, 1 = one or more failed.\n");
    printf("Running test_barrier_control\n");
    returnCode = ::main_test_barrier_control(argc, argv);
    if (returnCode != 0)
        return returnCode;

    printf("Running test_racing_lp\n");
    returnCode = ::main_test_racing_lp();
    if (returnCode != 0)
        return returnCode;

    return 0;
}

int main(int argc, char** argv)
{
    /*
    Suppress an popup window that Windows shows in response to a crash. See
    note at head of file.
    */
    WindowsErrorPopupBlocker();

    // Set code page to UTF-8 to avoid problems with non-ASCII characters in the output
    // #include <cstdlib>
    //system("chcp 65001 >NUL"); 

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--help" || arg == "-h") {
            printf(
                "ClpMakeTest - Runs tests for Clp aligned to 'make test'.\n"
                "\n"
                "Usage: ClpMakeTest [--help]\n"
                "\n"
                "Exit code: 0 = all checks passed, 1 = one or more checks failed.\n");
            return 0;
        }
        fprintf(stderr, "ClpMakeTest: unrecognized argument '%s' "
            "(try --help)\n",
            arg.c_str());
        return 1;
    }

    return test(argc, argv);
}