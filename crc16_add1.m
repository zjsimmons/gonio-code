function [OldCRC] = crc16_add1(data,OldCRC)

OldCRC = uint16(bitor(bitshift(OldCRC,-8,'uint16'),bitshift(OldCRC,8,'uint16')));
OldCRC = uint16(bitxor(OldCRC,data));
OldCRC = uint16(bitxor(OldCRC,bitshift(uint16(bitand(OldCRC,hex2dec('ff'))),-4,'uint16')));
OldCRC = uint16(bitxor(OldCRC,bitshift(bitshift(OldCRC,8,'uint16'),4,'uint16')));
OldCRC = uint16(bitxor(OldCRC,bitshift(bitshift(bitand(OldCRC,hex2dec('ff')),4,'uint16'),1,'uint16')));

end

