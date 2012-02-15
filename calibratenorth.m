 
north = zeros(8,3);
for i=1:8
    [imu, imup, file] = wsdread(strcat('WASP_INS_data/20120209/w315-20120209-',num2str(i,'%05d'),'.wsd'));
    north(i,:) = mean(imu(:,11:13));
end;

biases = zeros(4,3);
for i=1:4
     biases(i,1:2) = (north(i,1:2)+north(i+1,1:2))/2;
end

%bias = mean(biases,1)
for i=1:8
     north(i,:) = north(i,:) - bias315';
end

for i=1:4
     north(i*2,:) = north(i*2,:).*[-1 -1 1];
end

north = mean(north,1)