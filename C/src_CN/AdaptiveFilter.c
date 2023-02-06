/*
 * @file AdaptiveFilter.c
 *  
 * Adaptive Filter实现了一个归一化的最小均方自适应滤波器。
 * 为“期望信号”输入或“错误信号”输入提供了例程。
 *
 * Created on: Apr 14, 2014
 * Author: John Bang
 */

/******************************************************************************/
/* include block */
#include "AdaptiveFilter.h"

/******************************************************************************/
/** local definitions **/
static void AdaptWeights(AfData *pData);
static double Filter(double input, AfData *pData);
static double SquaredNorm(double *x, unsigned int length);

/******************************************************************************
 * AdaptiveFilterRun
 *
 * @param[in]     input  input signal sample
 * @param[in]     desired 信号样本期望
 *
 * @returns       adaptive filter output (估计信号期望)
 *
 * @note          运行归一化最小均方自适应滤波器并计算一个新的输出。
 *
 * @warning       none
 */
double AdaptiveFilterRun(double input, double desired, AfData *pData) {
	double output;

	output = Filter(input, pData);
	pData->Error = desired - output; 
	AdaptWeights(pData); 

	return output;
}



/******************************************************************************
 * AdaptiveFilterRunErrorIn
 *
 * @param[in]     input  input signal sample
 * @param[in]     error error signal sample (desired - output)
 * @param[in,out] pData  pointer to AdaptiveFilter parameter/state struct
 *
 * @returns       自适应滤波器输出(估计所需信号)
 *
 * @note          Runs the normalized least mean square adaptive filter and
 *  computes a new output.
 *
 * @warning       none
 */
double AdaptiveFilterRunErrorIn(double input, double error, AfData *pData) {
	double output;

	pData->Error = error; /* update the error */
	AdaptWeights(pData); /* update adaptive filter weights */
	output = Filter(input, pData); /* filter the input */

	return output;
}

/** internal functions **/

/***************************************************************************//**
* AdaptWeights
* 
* @param[in,out]     pData 指向自适应滤波器参数/状态结构体的指针
*
* @returns       none
* 
* @note          使用规范化最小均方算法更新pData->pWeights中的过滤器权重
* 
* @warning       none
*******************************************************************************/
static void AdaptWeights(AfData *pData) {
	double sn, normStepSize;
	int i;

	sn = SquaredNorm(pData->pBuffer,pData->Length); /* 计算规范项目 */
	normStepSize = (pData->StepSize)/(pData->Regularization + sn); /* 归一化步长 */

	for ( i = pData->Length - 1; i >= 0; i--) {
        /* wrap index */
        if (pData->BufferIdx >= pData->Length) {
            pData->BufferIdx = 0;
        }
        /* 归一化最小均方根方程 */
		pData->pWeights[i] += normStepSize * (pData->Error) * (pData->pBuffer[pData->BufferIdx++]);
	}
}


/***************************************************************************//**
* Filter
* 
* @param[in]     input 信号样本输入
* @param[in,out]     pData 自适应滤波器参数/状态结构体指针
*
* @returns       新的滤波器输出
* 
* @note          使用输入和当前滤波器的权重值计算出新的输出样本
* 
* @warning       none
*******************************************************************************/
static double Filter(double input, AfData *pData) {
	double output = 0;
	int i;
    
    /* 换行索引 */
    if (pData->BufferIdx >= pData->Length) {
        pData->BufferIdx = 0;
    }
    /* 利用新的输入覆盖旧的输入 */
	pData->pBuffer[pData->BufferIdx++] = input;

	for (i = pData->Length - 1; i >= 0; i--) {
        /* 换行索引 */
        if (pData->BufferIdx >= pData->Length) {
            pData->BufferIdx = 0;
        }
        /* 计算权重向量和缓冲区的内积 */
		output += (pData->pWeights[i]) * (pData->pBuffer[pData->BufferIdx++]);
	}
    
	return output;
}
/* End of Filter()*/
/******************************************************************************/

/***************************************************************************//**
* SquaredNorm
* 
* @param[in]     pInput 指向输入的指针
* @param[in]     length 缓冲区长度
*
* @returns       输入缓冲区的L2范数的平方
* 
* @note          计算输入缓冲区的L2范数平方，即每个元素的平方之和。
* 
* @warning       none
*******************************************************************************/
static double SquaredNorm(double *pInput, const unsigned int length) {
	double output = 0;
	unsigned int i;

    //注意:这可以在计算上得到改进，但目前这样做是为了明确地避免由于浮点加法的性质而造成的数值误差积累。
	for ( i = 0; i < length; i++ ) {
		output += pInput[i]*pInput[i]; 
	}
    
	return output;
}



