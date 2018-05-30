function [ As ] = tsmooth2( A, smoothlength )
%TSMOOTH Summary of this function goes here
%   Detailed explanation goes here
    [M,N,P]=size(A);
    As=A;
    for x=1:M
        As(x,:,:)=filt2d(reshape(A(x,:,:),N,P),'avgnans',1,smoothlength);
    end

end