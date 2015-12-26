function out = B_bar(B)
%The input B is the matrix containing the current location of all of the
%birds
global N
sumx = B(1,1);
sumy = B(2,1);
for n = 2:N
    sumx = sumx + B(1,N);
    sumy = sumy + B(2,N);
end
out = [sumx;sumy]*(1/N);