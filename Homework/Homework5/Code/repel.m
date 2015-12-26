function out = repel(B_k,B_l)
global rho delta N

difference = zeros(2,N-1);
mag_diff = zeros(1,N-1);
for l = 1:N-1
    difference(:,l) = B_k - B_l(:,l);
    xandy = difference(:,l);
    mag_diff(l) = sqrt(xandy(1)^2+xandy(2)^2);
end
sort_mag = sort(mag_diff,'descend');
sumof = zeros(2,N-1);
for close = 1:8
    locations = find(mag_diff == sort_mag(close));
    if size(B_k) == size(B_l(:,locations))
        sumof(:,close) = rho*((B_k - B_l(:,locations))./((B_k - B_l(:,locations)).^2+delta));
    else
        sumof(:,close) = [0;0];
    end
end
out = sum(sumof,2);
