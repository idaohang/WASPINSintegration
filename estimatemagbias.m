 
reading = zeros(6,3);
for i=1:6
    [imu, imup, file] = wsdread(strcat('WASP_INS_data/20120209/w220-20120209-',num2str(i,'%05d'),'.wsd'));
    reading(i,:) = mean(imu(:,11:13));
end;

reading
bias220 = [(reading(1,1) + reading(2,1))/2; (reading(3,2) + reading(4,2))/2; (reading(5,3) + reading(6,3))/2]

for i=1:6
    [imu, imup, file] = wsdread(strcat('WASP_INS_data/20120209/w315-20120209-',num2str(8+i,'%05d'),'.wsd'));
    reading(i,:) = mean(imu(:,11:13));
end;

reading
bias315 = [(reading(1,1) + reading(2,1))/2; (reading(3,2) + reading(4,2))/2; (reading(5,3) + reading(6,3))/2]