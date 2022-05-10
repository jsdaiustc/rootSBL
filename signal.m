function X=signal_sparse(M, alpha, SNR, K)
N_alpha=length(alpha);
A=exp(-i*pi*(0:M-1)'*sin(alpha*pi/180));
Vj=diag(sqrt((   10.^(SNR/10)   )/2));
S=Vj*(randn(N_alpha,K)+i*randn(N_alpha,K));
noise=sqrt(1/2)*(randn(M,K)+i*randn(M,K));
X=A*S+noise;
