#ifndef ADAPTIVEFILTER_H_
#define ADAPTIVEFILTER_H_

/* 自适应滤波器的参数结构体
 */
typedef struct {
	const double StepSize; /* 自适应滤波器步长 */
    const double Regularization; /* 正则化常数 */
	const unsigned int Length; /* 滤波器长度 */
	double *pBuffer; /* 指向输入缓冲区的指针 */
    unsigned int BufferIdx; /* 循环索引到输入缓冲区 */
	double *pWeights; /* 指向自适应滤波器权重的指针 */
	double Error; /* 指向输出错误（期望-输出）状态 */
} AfData;

double AdaptiveFilterRun (double input, double desired, AfData *pData);
double AdaptiveFilterRunErrorIn(double input, double error, AfData *pData);

#endif
