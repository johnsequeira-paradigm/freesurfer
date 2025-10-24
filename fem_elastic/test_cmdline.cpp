/**
 * Unit tests for command-line parsing in surf2vol
 * Tests the getopt-based argument parsing
 */

#include <iostream>
#include <cassert>
#include <cstring>
#include <vector>
#include <getopt.h>

// Simple test framework macros
#define TEST(name) void test_##name()
#define RUN_TEST(name) do { \
    std::cout << "Running " #name "..." << std::flush; \
    test_##name(); \
    std::cout << " PASSED" << std::endl; \
    tests_passed++; \
} while(0)

#define ASSERT_TRUE(expr) do { \
    if (!(expr)) { \
        std::cerr << "\n  ASSERTION FAILED: " #expr \
                  << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        exit(1); \
    } \
} while(0)

#define ASSERT_EQ(a, b) do { \
    if ((a) != (b)) { \
        std::cerr << "\n  ASSERTION FAILED: " #a " == " #b \
                  << "\n  Expected: " << (b) \
                  << "\n  Got: " << (a) \
                  << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        exit(1); \
    } \
} while(0)

#define ASSERT_STREQ(a, b) do { \
    if (std::string(a) != std::string(b)) { \
        std::cerr << "\n  ASSERTION FAILED: " #a " == " #b \
                  << "\n  Expected: " << (b) \
                  << "\n  Got: " << (a) \
                  << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        exit(1); \
    } \
} while(0)

int tests_passed = 0;

// Mock IoParams structure
struct IoParams {
    typedef std::vector<std::string> StringVectorType;

    StringVectorType vstrFixedSurf;
    StringVectorType vstrAparc;
    bool hasAparc;

    std::string strFixedMri;
    std::string strAseg;
    StringVectorType vstrMovingSurf;
    std::string strMovingMri;

    std::string strOutput;
    std::string strOutputField;
    std::string strOutputMesh;
    std::string strOutputSurf;
    std::string strOutputSurfAffine;
    std::string strGcam;
    std::string strOutputAffine;
    std::string dbgOutput;

    double eltVolMin, eltVolMax;
    float poissonRatio;
    float YoungModulus;
    bool compress;
    std::string strTransform;
    int iSteps;
    int iEndStep;
    double surfSubsample;
    std::string strDebug;
    bool bUseOldTopologySolver;
    bool bUsePialForSurf;

    IoParams() : hasAparc(false), eltVolMin(2), eltVolMax(21),
                 poissonRatio(0.3f), YoungModulus(10), compress(false),
                 iSteps(1), iEndStep(-1), surfSubsample(-1),
                 bUseOldTopologySolver(false), bUsePialForSurf(false),
                 strOutput("out.mgz"), strOutputField("out_field.mgz") {}

    int parse(int argc, char* argv[], std::string& errMsg);
};

// Simplified parse function for testing
int IoParams::parse(int argc, char* argv[], std::string& errMsg) {
    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"fixed-mri", required_argument, 0, 0},
        {"moving-mri", required_argument, 0, 0},
        {"aseg", required_argument, 0, 0},
        {"fixed-surf", required_argument, 0, 0},
        {"moving-surf", required_argument, 0, 0},
        {"aparc", required_argument, 0, 0},
        {"out", required_argument, 0, 'o'},
        {"out-field", required_argument, 0, 0},
        {"out-mesh", required_argument, 0, 0},
        {"elt-vol", required_argument, 0, 0},
        {"elt-vol-min", required_argument, 0, 0},
        {"elt-vol-max", required_argument, 0, 0},
        {"poisson", required_argument, 0, 0},
        {"young", required_argument, 0, 0},
        {"fem-steps", required_argument, 0, 0},
        {"compress", no_argument, 0, 0},
        {"topology-old", no_argument, 0, 0},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;

    // Reset getopt
    optind = 1;

    while ((c = getopt_long(argc, argv, "ho:", long_options, &option_index)) != -1) {
        if (c == 0) {
            std::string opt_name = long_options[option_index].name;

            if (opt_name == "fixed-mri")
                strFixedMri = optarg;
            else if (opt_name == "moving-mri")
                strMovingMri = optarg;
            else if (opt_name == "aseg")
                strAseg = optarg;
            else if (opt_name == "fixed-surf")
                vstrFixedSurf.push_back(optarg);
            else if (opt_name == "moving-surf")
                vstrMovingSurf.push_back(optarg);
            else if (opt_name == "aparc") {
                hasAparc = true;
                vstrAparc.push_back(optarg);
            }
            else if (opt_name == "out-field")
                strOutputField = optarg;
            else if (opt_name == "out-mesh")
                strOutputMesh = optarg;
            else if (opt_name == "elt-vol")
                eltVolMin = eltVolMax = atof(optarg);
            else if (opt_name == "elt-vol-min")
                eltVolMin = atof(optarg);
            else if (opt_name == "elt-vol-max")
                eltVolMax = atof(optarg);
            else if (opt_name == "poisson")
                poissonRatio = atof(optarg);
            else if (opt_name == "young")
                YoungModulus = atof(optarg);
            else if (opt_name == "fem-steps")
                iSteps = atoi(optarg);
            else if (opt_name == "compress")
                compress = true;
            else if (opt_name == "topology-old")
                bUseOldTopologySolver = true;
        }
        else if (c == 'h') {
            return -1;  // Help requested
        }
        else if (c == 'o') {
            strOutput = optarg;
        }
    }

    // Validation
    if (strFixedMri.empty())
        errMsg += "No fixed volume specified\n";
    if (strMovingMri.empty())
        errMsg += "No moving volume specified\n";
    if (vstrFixedSurf.empty())
        errMsg += "No fixed surface specified\n";
    if (vstrMovingSurf.empty())
        errMsg += "No moving surface specified\n";

    return errMsg.empty() ? 0 : 1;
}

// Helper to create argv
char** make_argv(const std::vector<std::string>& args) {
    char** argv = new char*[args.size()];
    for (size_t i = 0; i < args.size(); i++) {
        argv[i] = strdup(args[i].c_str());
    }
    return argv;
}

void free_argv(char** argv, int argc) {
    for (int i = 0; i < argc; i++) {
        free(argv[i]);
    }
    delete[] argv;
}

// Test cases

TEST(basic_required_args) {
    std::vector<std::string> args = {
        "surf2vol",
        "--fixed-mri", "fixed.mgz",
        "--moving-mri", "moving.mgz",
        "--fixed-surf", "lh.white",
        "--moving-surf", "lh.white.moved"
    };

    char** argv = make_argv(args);
    IoParams params;
    std::string errMsg;

    int result = params.parse(args.size(), argv, errMsg);

    ASSERT_EQ(result, 0);
    ASSERT_STREQ(params.strFixedMri.c_str(), "fixed.mgz");
    ASSERT_STREQ(params.strMovingMri.c_str(), "moving.mgz");
    ASSERT_EQ(params.vstrFixedSurf.size(), 1);
    ASSERT_STREQ(params.vstrFixedSurf[0].c_str(), "lh.white");
    ASSERT_EQ(params.vstrMovingSurf.size(), 1);
    ASSERT_STREQ(params.vstrMovingSurf[0].c_str(), "lh.white.moved");

    free_argv(argv, args.size());
}

TEST(missing_required_args) {
    std::vector<std::string> args = {
        "surf2vol",
        "--fixed-mri", "fixed.mgz"
    };

    char** argv = make_argv(args);
    IoParams params;
    std::string errMsg;

    int result = params.parse(args.size(), argv, errMsg);

    ASSERT_EQ(result, 1);  // Should fail
    ASSERT_TRUE(!errMsg.empty());
    ASSERT_TRUE(errMsg.find("moving volume") != std::string::npos);

    free_argv(argv, args.size());
}

TEST(optional_output_args) {
    std::vector<std::string> args = {
        "surf2vol",
        "--fixed-mri", "fixed.mgz",
        "--moving-mri", "moving.mgz",
        "--fixed-surf", "lh.white",
        "--moving-surf", "lh.white.moved",
        "-o", "custom_out.mgz",
        "--out-mesh", "mesh.tm3d"
    };

    char** argv = make_argv(args);
    IoParams params;
    std::string errMsg;

    int result = params.parse(args.size(), argv, errMsg);

    ASSERT_EQ(result, 0);
    ASSERT_STREQ(params.strOutput.c_str(), "custom_out.mgz");
    ASSERT_STREQ(params.strOutputMesh.c_str(), "mesh.tm3d");

    free_argv(argv, args.size());
}

TEST(fem_parameters) {
    std::vector<std::string> args = {
        "surf2vol",
        "--fixed-mri", "fixed.mgz",
        "--moving-mri", "moving.mgz",
        "--fixed-surf", "lh.white",
        "--moving-surf", "lh.white.moved",
        "--poisson", "0.35",
        "--young", "50",
        "--fem-steps", "5",
        "--elt-vol-min", "1.5",
        "--elt-vol-max", "25.0"
    };

    char** argv = make_argv(args);
    IoParams params;
    std::string errMsg;

    int result = params.parse(args.size(), argv, errMsg);

    ASSERT_EQ(result, 0);
    ASSERT_EQ(params.poissonRatio, 0.35f);
    ASSERT_EQ(params.YoungModulus, 50.0f);
    ASSERT_EQ(params.iSteps, 5);
    ASSERT_EQ(params.eltVolMin, 1.5);
    ASSERT_EQ(params.eltVolMax, 25.0);

    free_argv(argv, args.size());
}

TEST(boolean_flags) {
    std::vector<std::string> args = {
        "surf2vol",
        "--fixed-mri", "fixed.mgz",
        "--moving-mri", "moving.mgz",
        "--fixed-surf", "lh.white",
        "--moving-surf", "lh.white.moved",
        "--compress",
        "--topology-old"
    };

    char** argv = make_argv(args);
    IoParams params;
    std::string errMsg;

    int result = params.parse(args.size(), argv, errMsg);

    ASSERT_EQ(result, 0);
    ASSERT_TRUE(params.compress);
    ASSERT_TRUE(params.bUseOldTopologySolver);

    free_argv(argv, args.size());
}

TEST(multiple_surfaces) {
    std::vector<std::string> args = {
        "surf2vol",
        "--fixed-mri", "fixed.mgz",
        "--moving-mri", "moving.mgz",
        "--fixed-surf", "lh.white",
        "--fixed-surf", "lh.pial",
        "--moving-surf", "lh.white.moved",
        "--moving-surf", "lh.pial.moved"
    };

    char** argv = make_argv(args);
    IoParams params;
    std::string errMsg;

    int result = params.parse(args.size(), argv, errMsg);

    ASSERT_EQ(result, 0);
    ASSERT_EQ(params.vstrFixedSurf.size(), 2);
    ASSERT_EQ(params.vstrMovingSurf.size(), 2);
    ASSERT_STREQ(params.vstrFixedSurf[0].c_str(), "lh.white");
    ASSERT_STREQ(params.vstrFixedSurf[1].c_str(), "lh.pial");

    free_argv(argv, args.size());
}

TEST(elt_vol_single_value) {
    std::vector<std::string> args = {
        "surf2vol",
        "--fixed-mri", "fixed.mgz",
        "--moving-mri", "moving.mgz",
        "--fixed-surf", "lh.white",
        "--moving-surf", "lh.white.moved",
        "--elt-vol", "10.0"
    };

    char** argv = make_argv(args);
    IoParams params;
    std::string errMsg;

    int result = params.parse(args.size(), argv, errMsg);

    ASSERT_EQ(result, 0);
    ASSERT_EQ(params.eltVolMin, 10.0);
    ASSERT_EQ(params.eltVolMax, 10.0);

    free_argv(argv, args.size());
}

TEST(default_values) {
    std::vector<std::string> args = {
        "surf2vol",
        "--fixed-mri", "fixed.mgz",
        "--moving-mri", "moving.mgz",
        "--fixed-surf", "lh.white",
        "--moving-surf", "lh.white.moved"
    };

    char** argv = make_argv(args);
    IoParams params;
    std::string errMsg;

    int result = params.parse(args.size(), argv, errMsg);

    ASSERT_EQ(result, 0);
    ASSERT_EQ(params.poissonRatio, 0.3f);
    ASSERT_EQ(params.YoungModulus, 10.0f);
    ASSERT_EQ(params.iSteps, 1);
    ASSERT_EQ(params.eltVolMin, 2.0);
    ASSERT_EQ(params.eltVolMax, 21.0);
    ASSERT_STREQ(params.strOutput.c_str(), "out.mgz");
    ASSERT_FALSE(params.compress);

    free_argv(argv, args.size());
}

// Main test runner
int main() {
    std::cout << "\n=== Running Command-Line Parsing Tests ===\n" << std::endl;

    RUN_TEST(basic_required_args);
    RUN_TEST(missing_required_args);
    RUN_TEST(optional_output_args);
    RUN_TEST(fem_parameters);
    RUN_TEST(boolean_flags);
    RUN_TEST(multiple_surfaces);
    RUN_TEST(elt_vol_single_value);
    RUN_TEST(default_values);

    std::cout << "\n=== All " << tests_passed << " tests passed! ===\n" << std::endl;
    return 0;
}
