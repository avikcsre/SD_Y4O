{{ USAGE OF SD_Y4O }}

[ SD_Y4O_1.m ]

** Input: **
  C3 matrix  in PolSARPro format (C11, C22, C33, ReC12, ImC12, C12, C21, ReC13, ImC13, C13, C31, ReC23, ImC23, C23, C32)
	Number of Rows
	Number of Cols
	Number of Looks (put arbitrary, e.g: 9)
	Enter Window size
	Angle increment Step
	
** Output: **
  THEETA_Real_Rot.bin, THEETA_HD.bin, HD_T33_T22.bin

[ SD_Y4O_2.m ]

** Input: **
  T3 matrix in PolSARPro format (T11, T22, T33, T12_real, T12_imag, T13_real, T13_imag, T23_real, T23_imag), 
  THEETA_HD.bin (from previous output)
  Enter Look Search Range
  Enter step size for L
  Output:  Max_Look.bin, HD_maxL.bin

___________________

alpha = (0.5/45)*| THEETA_HD.bin | + 0.5
beta = 1 – alpha

The outputs of the powers (P_s, P_d, P_v, P_c) are obtained from PolSARPro for Y4O Decomposition 
The Final Outputs:
P_v^n = P_v*(1 - HD_maxL.bin)
P_d^n = P_d + alpha*P_v* HD_maxL.bin
P_s^n = P_s + beta*P_v* HD_maxL.bin
P_c^n = P_c


