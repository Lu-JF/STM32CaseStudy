#ifndef __HILBERT_H
#define __HILBERT_H

void abs_hilbert(int n,double x[]);
void hilbert(int n,double x[]);
void fht(int n,double x[]);
void get_hilbert_angle(int n,double x[]);
void unwrap(int n,double x[]);
void resample(int n,double x[],double y[],double time[]);
void dfft(int n,int k,double fr[],double fi[],double pr[],double pi[]);






#endif





