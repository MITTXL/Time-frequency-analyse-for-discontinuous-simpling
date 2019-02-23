时间：20190106

实验问题：对于均匀采样信号的FFT和LS验证。

出现的问题：按照论文中的思路进行编程时，出现LS方法和FFT效果相差甚远的结果。对比学长实验程序，其原因在于变换基的产生的顺序，现详细记录如下。

变换基产生程序为：

![image1](C:\Users\Aska Hu\OneDrive\Time-Frequencyrate\2019.01 非均匀时频变换\FFT+LS\1_FFT&LS实验复现\image1.png)

在我们实验程序中，直流分量在矩阵的第一列

![image2](C:\Users\Aska Hu\OneDrive\Time-Frequencyrate\2019.01 非均匀时频变换\FFT+LS\1_FFT&LS实验复现\image2.png)

更改后将 直流分量放置在最后一列：结果正确

![image3](C:\Users\Aska Hu\OneDrive\Time-Frequencyrate\2019.01 非均匀时频变换\FFT+LS\1_FFT&LS实验复现\image3.png)

将直流分量放在最后一行解决问题

![image5](C:\Users\Aska Hu\OneDrive\Time-Frequencyrate\2019.01 非均匀时频变换\FFT+LS\1_FFT&LS实验复现\image5.png)

运行结果如下![image6](C:\Users\Aska Hu\OneDrive\Time-Frequencyrate\2019.01 非均匀时频变换\FFT+LS\1_FFT&LS实验复现\image6.png)