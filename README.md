# CPU Scheduling Simulator

## Overview
This project is a CPU scheduling simulator implemented in C++ as part of the "Operating Systems" course (CC373) at Alexandria University. The program analyzes and visualizes various CPU scheduling policies, providing statistical output and execution traces.

## Supported Scheduling Algorithms
The following CPU scheduling algorithms are implemented:
1. **FCFS** (First Come First Serve)
2. **RR** (Round Robin)
3. **SPN** (Shortest Process Next)
4. **SRT** (Shortest Remaining Time)
5. **HRRN** (Highest Response Ratio Next)
6. **FB-1** (Feedback with q=1 for all queues)
7. **FB-2i** (Feedback with q=2^i for each queue)
8. **Aging**

## Installation and Compilation
To compile the project, use the provided `Makefile`. Run the following command in the project directory:

```sh
make
```

This will generate an executable named `lab6` in the same directory.

## Usage
The program reads input from `stdin` and writes output to `stdout`. Input can be redirected from a file or piped:

```sh
./lab6 < input.txt
```

or

```sh
cat input.txt | ./lab6
```

### Input Format
The input consists of multiple lines:
1. **Line 1**: `trace` or `stats`
   - `trace`: Visualizes process execution over time.
   - `stats`: Outputs scheduling statistics.
2. **Line 2**: Comma-separated list of scheduling policies to analyze. Some policies require parameters (e.g., `2-4` for Round Robin with quantum = 4).
3. **Line 3**: Simulation end time (integer).
4. **Line 4**: Number of processes (integer).
5. **Lines 5+**: Process details in comma-separated format:
   - For algorithms **1-7**: `Process Name, Arrival Time, Service Time`
   - For **Aging (8)**: `Process Name, Arrival Time, Priority`

### Example Input
```
trace
1,2-4,3
50
3
A,0,5
B,2,3
C,4,2
```

### Output Format
- **Trace Mode**: Uses `*` for running processes and `.` for ready processes.
- **Stats Mode**: Displays performance metrics for scheduled processes.

To verify correctness, compare output with expected test cases:
```sh
cat test_output.txt | diff expected_output.txt -
```
---
**Author**: Marly Magued 
**Course**: Operating Systems (CC373)  
**Institution**: Alexandria University, Faculty of Engineering  
