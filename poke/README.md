# POKÉ with Dual HD isogenies

## Our Contributions
This project is built upon the original C-implementation of POKÉ, INKE, and PIKE (available at [https://github.com/gusgkr0117/poke](https://github.com/gusgkr0117/poke); the original README is also appended below). In this artifact, we introduce and implement five new protocol variants as proposed in our paper.

To facilitate benchmarking and evaluation, we have implemented the **NIST Level 1** security parameters for all our new variants. The mapping between the source code directories and the schemes described in our paper is as follows:

| Source Code Directory | Protocol Variant in Paper |
| :--- | :--- |
| `poket/` | **POKÉD** |
| `piket/` | **PIKED** |
| `pokeb_p/` | **B-POKÉD** |
| `pokeb_m/` | **B-POKÉD** (alternative implementation) |
| `pikeb/` | **B-PIKED** |

## Build and Test

To build the project and run the performance benchmarks, ensure that `Clang` and `GMP` are installed. The following commands will build the project and execute the Level 1 tests for both the original schemes and our new variants, allowing for direct efficiency comparisons.

```bash
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_BUILD_TYPE=Release ..
make
ctest -V -R "poke-test*" -E "hard|lvl3|lvl5"
```
*(Note: The `-E "hard|lvl3|lvl5"` flag excludes the Level 3, Level 5, and "hard" tests, isolating the Level 1 implementations for a focused comparison.)*

## Performance Benchmark Example

The following is a sample test output run on an **Intel Xeon w3-2423** CPU. To ensure accurate and stable benchmarking, CPU features such as Turbo Boost and Hyper-Threading were disabled during the evaluation.

```text
Test project ./build
Constructing a list of tests
Done constructing a list of tests
Updating test list for fixtures
Added 0 tests to meet fixture requirements
Checking test dependency graph...
Checking test dependency graph end
test 13
    Start 13: poke-test_poke_lvl1

13: Test command: ./build/src/poke/ref/poke_lvl1/test/poke-test_poke_lvl1 "100"
13: Working Directory: ./build/src/poke/ref/poke_lvl1/test
13: Test timeout computed to be: 1500
13: test loops : 100
13:   keygen takes .................................... 244439311.320000 cycles
13:   encrypt takes .................................... 56558466.020000 cycles
13:   decrypt takes .................................... 82454150.680000 cycles
1/7 Test #13: poke-test_poke_lvl1 ..............   Passed   18.17 sec
test 31
    Start 31: poke-test_pike_lvl1

31: Test command: ./build/src/poke/ref/pike_lvl1/test/poke-test_pike_lvl1 "100"
31: Working Directory: ./build/src/poke/ref/pike_lvl1/test
31: Test timeout computed to be: 1500
31: test loops : 100
31:   keygen takes .................................... 261885914.120000 cycles
31:   encrypt takes .................................... 50365142.160000 cycles
31:   decrypt takes .................................... 45384336.800000 cycles
2/7 Test #31: poke-test_pike_lvl1 ..............   Passed   16.95 sec
test 37
    Start 37: poke-test_pikeb_lvl1

37: Test command: ./build/src/poke/ref/pikeb_lvl1/test/poke-test_pikeb_lvl1 "100"
37: Working Directory: ./build/src/poke/ref/pikeb_lvl1/test
37: Test timeout computed to be: 1500
37: test loops : 100
37:   keygen takes .................................... 75882657.620000 cycles
37:   encrypt takes .................................... 99729725.360000 cycles
37:   decrypt takes .................................... 13231732.840000 cycles
3/7 Test #37: poke-test_pikeb_lvl1 .............   Passed    8.95 sec
test 39
    Start 39: poke-test_pokeb_p_lvl1

39: Test command: ./build/src/poke/ref/pokeb_p_lvl1/test/poke-test_pokeb_p_lvl1 "100"
39: Working Directory: ./build/src/poke/ref/pokeb_p_lvl1/test
39: Test timeout computed to be: 1500
39: test loops : 100
39:   keygen takes .................................... 79377759.200000 cycles
39:   encrypt takes .................................... 97824254.460000 cycles
39:   decrypt takes .................................... 12481427.440000 cycles
4/7 Test #39: poke-test_pokeb_p_lvl1 ...........   Passed    8.99 sec
test 42
    Start 42: poke-test_pokeb_m_lvl1

42: Test command: ./build/src/poke/ref/pokeb_m_lvl1/test/poke-test_pokeb_m_lvl1 "100"
42: Working Directory: ./build/src/poke/ref/pokeb_m_lvl1/test
42: Test timeout computed to be: 1500
42: test loops : 100
42:   keygen takes .................................... 82351950.180000 cycles
42:   encrypt takes .................................... 99409835.760000 cycles
42:   decrypt takes .................................... 12545388.020000 cycles
5/7 Test #42: poke-test_pokeb_m_lvl1 ...........   Passed    9.21 sec
test 45
    Start 45: poke-test_poket_lvl1

45: Test command: ./build/src/poke/ref/poket_lvl1/test/poke-test_poket_lvl1 "100"
45: Working Directory: ./build/src/poke/ref/poket_lvl1/test
45: Test timeout computed to be: 1500
45: test loops : 100
45:   keygen takes .................................... 61562530.960000 cycles
45:   encrypt takes .................................... 26977589.880000 cycles
45:   decrypt takes .................................... 23271424.820000 cycles
6/7 Test #45: poke-test_poket_lvl1 .............   Passed    5.30 sec
test 48
    Start 48: poke-test_piket_lvl1

48: Test command: ./build/src/poke/ref/piket_lvl1/test/poke-test_piket_lvl1 "100"
48: Working Directory: ./build/src/poke/ref/piket_lvl1/test
48: Test timeout computed to be: 1500
48: test loops : 100
48:   keygen takes .................................... 62739929.940000 cycles
48:   encrypt takes .................................... 27952907.980000 cycles
48:   decrypt takes .................................... 25540954.320000 cycles
7/7 Test #48: poke-test_piket_lvl1 .............   Passed    5.51 sec

The following tests passed:
    poke-test_poke_lvl1
    poke-test_pike_lvl1
    poke-test_pikeb_lvl1
    poke-test_pokeb_p_lvl1
    poke-test_pokeb_m_lvl1
    poke-test_poket_lvl1
    poke-test_piket_lvl1

100% tests passed, 0 tests failed out of 7

Total Test time (real) =  73.09 sec
```

---
# Original Project README Below

# 🥗POKÉ : POint-based Key Exchange and ✒️INKE : INtermediate-curve-based Key Exchange
This implementation is based on [SQIsign](https://github.com/SQISign/sqisign2d-west-ac24) (Apache License 2.0) and has been modified.

(2026.03.13 update) The [PIKE](https://eprint.iacr.org/2026/473) code is optimized and based on [this repository](https://github.com/Kaizhan-Lin/PIKE-C-Implementation).
+ C-Implementation of the most efficient isogeny-based PKE protocols
+ GMP must be installed
+ Clang should be used

# How-to-use
## Install Clang
```bash
sudo bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"
```

## Run POKE
```bash
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_BUILD_TYPE=Release ..
make precomp
make
ctest -V -R "poke-test*" -E "hard"
```
## Output
Run in Apple M2 CPU 3.49GHz
```
test 13
    Start 13: poke-test_poke_lvl1

13: Test command: poke/ref/poke_lvl1/test/poke-test_poke_lvl1 "100"
13: Working Directory: poke/build/src/poke/ref/poke_lvl1/test
13: Test timeout computed to be: 1500
13: test loops : 100
13:   keygen takes .................................... 53593587.200000 nsec
13:   encrypt takes .................................... 12392670.720000 nsec
13:   decrypt takes .................................... 17681781.760000 nsec
1/9 Test #13: poke-test_poke_lvl1 ..............   Passed    8.85 sec
test 16
    Start 16: poke-test_poke_lvl3

16: Test command: poke/ref/poke_lvl3/test/poke-test_poke_lvl3 "100"
16: Working Directory: poke/build/src/poke/ref/poke_lvl3/test
16: Test timeout computed to be: 1500
16: test loops : 100
16:   keygen takes .................................... 191975690.240000 nsec
16:   encrypt takes .................................... 44050311.680000 nsec
16:   decrypt takes .................................... 64517629.440000 nsec
2/9 Test #16: poke-test_poke_lvl3 ..............   Passed   30.30 sec
test 19
    Start 19: poke-test_poke_lvl5

19: Test command: poke/ref/poke_lvl5/test/poke-test_poke_lvl5 "100"
19: Working Directory: poke/build/src/poke/ref/poke_lvl5/test
19: Test timeout computed to be: 1500
19: test loops : 100
19:   keygen takes .................................... 392021795.840000 nsec
19:   encrypt takes .................................... 89838517.760000 nsec
19:   decrypt takes .................................... 130401054.720000 nsec
3/9 Test #19: poke-test_poke_lvl5 ..............   Passed   61.49 sec
test 22
    Start 22: poke-test_inke_lvl1

22: Test command: poke/ref/inke_lvl1/test/poke-test_inke_lvl1 "100"
22: Working Directory: poke/build/src/poke/ref/inke_lvl1/test
22: Test timeout computed to be: 1500
22: test loops : 100
22:   keygen takes .................................... 44252917.760000 nsec
22:   encrypt takes .................................... 11813091.840000 nsec
22:   decrypt takes .................................... 8816678.400000 nsec
4/9 Test #22: poke-test_inke_lvl1 ..............   Passed    6.86 sec
test 25
    Start 25: poke-test_inke_lvl3

25: Test command: poke/ref/inke_lvl3/test/poke-test_inke_lvl3 "100"
25: Working Directory: poke/build/src/poke/ref/inke_lvl3/test
25: Test timeout computed to be: 1500
25: test loops : 100
25:   keygen takes .................................... 131215429.120000 nsec
25:   encrypt takes .................................... 36028861.440000 nsec
25:   decrypt takes .................................... 26732697.600000 nsec
5/9 Test #25: poke-test_inke_lvl3 ..............   Passed   19.63 sec
test 28
    Start 28: poke-test_inke_lvl5

28: Test command: poke/ref/inke_lvl5/test/poke-test_inke_lvl5 "100"
28: Working Directory: poke/build/src/poke/ref/inke_lvl5/test
28: Test timeout computed to be: 1500
28: test loops : 100
28:   keygen takes .................................... 281477163.520000 nsec
28:   encrypt takes .................................... 76106129.920000 nsec
28:   decrypt takes .................................... 56573406.720000 nsec
6/9 Test #28: poke-test_inke_lvl5 ..............   Passed   41.65 sec
test 31
    Start 31: poke-test_pike_lvl1

31: Test command: poke/ref/pike_lvl1/test/poke-test_pike_lvl1 "100"
31: Working Directory: poke/build/src/poke/ref/pike_lvl1/test
31: Test timeout computed to be: 1500
31: test loops : 100
31:   keygen takes .................................... 57891683.840000 nsec
31:   encrypt takes .................................... 11052577.280000 nsec
31:   decrypt takes .................................... 9862195.200000 nsec
7/9 Test #31: poke-test_pike_lvl1 ..............   Passed    8.25 sec
test 33
    Start 33: poke-test_pike_lvl3

33: Test command: poke/ref/pike_lvl3/test/poke-test_pike_lvl3 "100"
33: Working Directory: poke/build/src/poke/ref/pike_lvl3/test
33: Test timeout computed to be: 1500
33: test loops : 100
33:   keygen takes .................................... 177139399.680000 nsec
33:   encrypt takes .................................... 34158737.920000 nsec
33:   decrypt takes .................................... 30239623.680000 nsec
8/9 Test #33: poke-test_pike_lvl3 ..............   Passed   24.40 sec
test 35
    Start 35: poke-test_pike_lvl5

35: Test command: poke/ref/pike_lvl5/test/poke-test_pike_lvl5 "100"
35: Working Directory: poke/build/src/poke/ref/pike_lvl5/test
35: Test timeout computed to be: 1500
35: test loops : 100
35:   keygen takes .................................... 376436930.560000 nsec
35:   encrypt takes .................................... 72581629.440000 nsec
35:   decrypt takes .................................... 63933672.960000 nsec
9/9 Test #35: poke-test_pike_lvl5 ..............   Passed   51.57 sec

The following tests passed:
	poke-test_poke_lvl1
	poke-test_poke_lvl3
	poke-test_poke_lvl5
	poke-test_inke_lvl1
	poke-test_inke_lvl3
	poke-test_inke_lvl5
	poke-test_pike_lvl1
	poke-test_pike_lvl3
	poke-test_pike_lvl5

100% tests passed, 0 tests failed out of 9

Total Test time (real) = 253.03 sec
```

# Reference
+ [INKE paper](https://eprint.iacr.org/2025/1458)
+ [PIKE paper](https://eprint.iacr.org/2026/473) and [C-implementation](https://github.com/Kaizhan-Lin/PIKE-C-Implementation)
+ [POKÉ paper](https://eprint.iacr.org/2024/624)
+ [SQISign git](https://github.com/SQISign/sqisign2d-west-ac24)
+ [fiat-crypto](https://github.com/mit-plv/fiat-crypto)