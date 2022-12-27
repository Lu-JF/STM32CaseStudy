#ifndef __FILTFILT_H
#define __FILTFILT_H


void iir_filt(double *x,int xlen,double *a,double *b,int nfilt,double *y);
void filtfilt(double *x,int xlen,double *a,double *b,int nfilt,double *y);






#endif


