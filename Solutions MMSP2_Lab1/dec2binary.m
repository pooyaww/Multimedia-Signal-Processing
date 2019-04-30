function out = dec2binary(y, Nbits)
% dec2binary converts decimal numbers into binary numbers
% y = dec2binary(x,n) converts the values of x into binary nubers on
% n bits

% First we create the LookUp Table for the 
% interger value <-> binary representation binding

for i = 1:2^Nbits       % cycle over all integer values between 0:2^Nbits
    int_val = i - 1;
    for j = Nbits:-1:1  % conversion of the actual value into Nbits bits
        lut(i, j) = mod(int_val, 2);    % remainder
        int_val = floor(int_val/2);     % rounded quotient
    end
end

% Finally we assign to each sample of y its binary representation
out = lut(y(:)+ 1, :);