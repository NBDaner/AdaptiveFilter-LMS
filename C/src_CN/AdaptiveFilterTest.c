 /*
 * @file AdaptiveFilterTest.c
 *  
* AdaptiveFilter测试程序:
* 1.创建自适应过滤器数据结构
* 2.创建一个固定的测试过滤器
* 3.生成一个随机输入信号
* 4.运行自适应过滤器以识别固定的测试过滤器权重
* 5.计算不对中和平方误差度量并打印到标准输出
* 6.根据预期的收敛阈值，报告通过/失败的标准输出
 *
 * Created on: Apr 14, 2014
 * Author: John Bang
 */

/******************************************************************************/
/* include block */
#include "AdaptiveFilter.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/** local definitions **/
static void InitWeights();
static double Filter(double input);
static double ComputeMisalignment();
static void PrintIterationStatus(unsigned int iteration);
static void PrintPassFailStatus();

/* Adaptive Filter parameter/state information ********************************/

/* 主要的测试参数 */
#define STEPSIZE (0.3) /* 自适应滤波器步长 */
#define REGULARIZATION (1.0E-10) /* 自适应滤波器正则化常数 */
#define NUM_TAPS (30) /* 自适应滤波器抽头数量 */
#define ITERATIONS (5000) /* 迭代次数 */
#define MISALIGNMENT_PASS_THRESH (-290.0) /* dB threshold for pass/fail test */
#define SQUARED_ERROR_PASS_THRESH (-290.0) /* dB threshold for pass/fail test */
#define DB_EPSILON (1.0E-40) /* allows minimum 10*log10() value of -400dB */
#define RAND_SEED (824) /* explicit random seed for test repeatability */

/* 测试状态 */
static double testWeights[NUM_TAPS];
static double testBuffer[NUM_TAPS];
static unsigned int testBufferIdx = 0;
static double squaredErrorDb, misalignmentDb;

/* 自适应滤波器数据 */
static double inBuffer[NUM_TAPS] = { 0 };
static double weights[NUM_TAPS] = { 0 };
static AfData Adata = {
		STEPSIZE,
        REGULARIZATION,
		NUM_TAPS,
		inBuffer,
        0,  //初始化索引
		weights,
		0.0 //初始化误差
};


/**
 * @param[in]     none
 *
 * @returns       none
 *
 * @note          在具有固定测试过滤器的系统中运行自适应过滤器，
 *                并根据上面列出的参数中定义的期望跟踪性能指标(错位和平方误差)。
 *
 * @warning       none
 */
void AdaptiveFilterTestRun() {
	double input, desired, output;
	unsigned int i;
    
    srand(RAND_SEED); /* set random seed for repeatability */

	InitWeights(); /* 初始化固定的测试滤波器 */

	for ( i = 0; i < ITERATIONS; i++) {
        /* 在区间(-1,1)上生成一个随机输入样本 */
		input = ( 2 * (double)rand() / (double)RAND_MAX ) - 1;
		desired = Filter(input); /* 运行固定的测试过滤器 */
		output = AdaptiveFilterRun(input, desired, &Adata);
        
        /* 计算性能指标 */
		squaredErrorDb = 10 * log10( DB_EPSILON + (Adata.Error) * (Adata.Error) );
        misalignmentDb = 10 * log10( DB_EPSILON + ComputeMisalignment() );
        
        PrintIterationStatus(i+1); /* print performance for this iteration */
	}
    PrintPassFailStatus(); /* print whether expected performance was acheived */

}

/** internal functions **/

/**
* @param[in]     none
*
* @returns       none
*
* @note          使用随机数初始化固定的测试过滤器权重
* 
* @warning       none
*******************************************************************************/
static void InitWeights() {
	unsigned int i;
	for ( i = 0; i < NUM_TAPS; i++) {
        /* 在区间(-1,1)上使用随机数初始化 */
		testWeights[i] = ( 2 * (double)rand() / (double)RAND_MAX ) - 1;
	}
}


/**
* @param[in]     input input signal sample
*
* @returns       new test filter output
* 
* @note          Computes a new output sample using the input and the fixed
*  test filter
* 
* @warning       none
*******************************************************************************/
static double Filter(double input) {
	double output = 0;
	int i;

    if (testBufferIdx >= NUM_TAPS) {
        testBufferIdx = 0;
    }
    testBuffer[testBufferIdx++] = input;
    
	for ( i = NUM_TAPS - 1; i >= 0; i--) {
        if (testBufferIdx >= NUM_TAPS) {
            testBufferIdx = 0;
        }
		output += testWeights[i] * testBuffer[testBufferIdx++];
	}

	return output;
}


/**
* @param[in]     none
*
* @returns       滤波器权数偏离
* 
* @note         计算测试滤波器和自适应滤波器之间的滤波器权重不对中，
                由测试滤波器的l2范数的平方归一化:(||Wtest - Wadaptive||^2) / (||Wtest||^2)
* 
* @warning       none
*******************************************************************************/
static double ComputeMisalignment() {
    unsigned int i;
    double difference;
    double diffSqrdNorm = 0.0, testSqrdNorm = 0.0;
    
    for ( i = 0; i < NUM_TAPS; i++) {
        difference = testWeights[i] - Adata.pWeights[i]; /* 权重差 */
        
        /* 累计平方项 */
        diffSqrdNorm += difference * difference;
        testSqrdNorm += testWeights[i] * testWeights[i];
    }
    
    return ( diffSqrdNorm / testSqrdNorm ); /* 返回归一化误差 */
}

/**
* @param[in]     iteration current iteration number
*
* @returns       none
* 
* @note          打印迭代相关的性能指标
* 
* @warning       none
*******************************************************************************/
static void PrintIterationStatus(unsigned int iteration) {
    printf("Iteration: %d\n", iteration);   // 迭代次数
    printf("Misalignment (dB): %f\n",misalignmentDb);  //偏差
    printf("Squared error (dB): %f\n",squaredErrorDb);  //均方误差
}

/**
* @param[in]     none
*
* @returns       none
* 
* @note          打印自适应滤波器系统的pass/fail状态
* 
* @warning       none
*******************************************************************************/
static void PrintPassFailStatus() {
    if (misalignmentDb > MISALIGNMENT_PASS_THRESH) {
        printf("FAIL: Misalignment !< %.0f\n",MISALIGNMENT_PASS_THRESH);
    }
    else {
        printf("PASS: Misalignment < %.0f\n",MISALIGNMENT_PASS_THRESH);
    }
    if (squaredErrorDb > SQUARED_ERROR_PASS_THRESH) {
        printf("FAIL: Squared Error !< %.0f\n",SQUARED_ERROR_PASS_THRESH);
    }
    else {
        printf("PASS: Squared Error < %.0f\n",SQUARED_ERROR_PASS_THRESH);
    }
}