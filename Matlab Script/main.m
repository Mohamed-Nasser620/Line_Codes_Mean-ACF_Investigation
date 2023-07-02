%Initialization of variables
A = 4;  %trasmitted voltage level
ensemble_size = 20000;  %number of realizations in the ensemble
num_bits = 100; %number o2f bits in each realization
bit_time = 70;  %pulse width to transmite one bit
sampling_time = 10; %DAC sampling time
num_samples_per_bit = floor(bit_time / sampling_time);  %numbers of samples in on bit (7)
num_samples = num_samples_per_bit*num_bits; %total number of samples in one realization

%Formation of the signals according to the required line code:
%fistly, generating the data as a matrix with nubmer of rows equal to the
%number of realization and number of coloumns represent the bits for each one adding a bit for delay 
data_bits = randi([0,1], ensemble_size, num_bits+1); 
%secondly, generate an empty matrix to store the data after transforming it
%form bits to samples
data_expanded = zeros(ensemble_size, num_samples+num_samples_per_bit);
%thirdly, Determining from the user the type of line coding to sample the bits
Line_code = input("1. Unipolar NRZ\n2. Polar NRZ\n3. Polar RZ\nEnter The Required Line Code: ");
%According to the input, the line code will be chosen
switch Line_code
    case 1  %Unipolar NRZ
        %transmit '0' as 0 and '1' as A
        data_symbols = data_bits * A;
        %transform bits to samples in all the realizations
        for i = 1:ensemble_size
            data_expanded(i,:) = repelem(data_symbols(i,:),num_samples_per_bit);
        end
    case 2  %Polar NRZ
        %transmit '0' as -A and '1' as A
        data_symbols = (2*data_bits-1)*A;
        %transform bits to samples in all the realizations
        for i = 1:ensemble_size
            data_expanded(i,:) = repelem(data_symbols(i,:),num_samples_per_bit);
        end
    case 3  %Polar RZ
        %transmit '0' as -A and '1' as A for the ceil of half the bit
        %period (4 out of 7 in our case)
        for i = 1:ensemble_size
            for j = 1:num_bits+1
                % Index of first sample for current bit
                start_idx = ((j - 1) * num_samples_per_bit) + 1;
                if data_bits(i,j) == 1 % Transmit '1' as +A followed by zero for the remaining bit period
                    data_expanded(i,start_idx:start_idx+floor(num_samples_per_bit/2)) = A;
                else % Transmit '0' as -A followed by zero for the remaining bit period
                    data_expanded(i,start_idx:start_idx+floor(num_samples_per_bit/2)) = -A;
                end
            end
        end       
end

%Adding a random dalay at the start of each realization
%firstly, generating a delay less than the number of samples per bit for all realizations
Td = randi([0,(num_samples_per_bit-1)],ensemble_size,1);
%defining a matrix to store the data after the delay
data = zeros(ensemble_size, num_samples);
%looping on all the realization and store the data from the end of the
%delay until getting all the samples
for i = 1:ensemble_size
    data(i,:)= data_expanded(i, Td(i)+1 : num_samples + Td(i));
end

%Ploting the final Realizations generated
figure('Name','Realizations');
%define a variable (t) represent the x-axis start from 0 and has a value with
%each sampling time until reaching the number of samples-1 as we start from 0 
t = 0:sampling_time:(num_samples-1)*sampling_time;
%looping 5 times to plot first 5 realization
for i = 1 : 5
    %plot the 5 realization in 1 figure inorder in a vertical position
    subplot(5,1,i);
    %set the value of x-axis as (t) and the corresponding y-axis is the generated data
    plot(t, data(i,:)); 
    %set a title for the graph
    str = sprintf("Realization %d",i);
    title(str);
    %set a limets from -6 to 6 on the y-axis
    ylim([-6 6]);
    %setting labels for x-axis and y-axis
    xlabel("Time(ms)");
    ylabel("Volts(V)");
end


%Q1: Getting the statistical Mean
%getting the sum of elements in each column and divided by their numeber which in the
%number of realization to get the mean
stat_mean = sum(data, 1) / ensemble_size;
%plotting the statistical mean
figure ('Name', 'Mean');
plot(t, stat_mean);
title ("Statistical Mean");
ylim([-6,6]);
xlabel("Time(ms)")
ylabel("Mean(V)");

%Q2: Determine if the random process is stationary
%Firstly, to say that the process is stationary, the mean must be constant
%which is previously checked
%Secondly, we will compare autocorrelations to check if it is dependent on 
%the time difference, so we have to get the autocorrelation in Q3 first  
%Q3: Statistical Autocorrelation calculations
%generate an empty matrix to store the result in it
stat_acf = zeros(size(data(1,1:end)));
%loop on the coulmns of the data matrix to multiply the first coloumn by all 
%the other coulmns and take their sum and divided by the number of elemnts
%in that coloumn which is represented by the ensemble size to get the right
%sided autocorrelation from 0 to the number of samples
for i = 1 : num_samples
    stat_acf(1,i) = sum((data(:,1) .* data(:,i))) / ensemble_size;
end
%get the left side by fliping the previously calculated autocorrelation and
%concatenate them to get the final statistical autocorrelation
stat_acf = cat (2, fliplr(stat_acf(2:num_samples)), stat_acf);

%Statistical Autocorrelation for another sample to check stationary
%choose a random number from the first quarter of our samples
random_check = randi([1,floor(num_samples/4)]);
%make an empty matrix to store the result
stat_acf2 = zeros(size(data(1,1:num_samples-random_check)));
%multiply from the random column with its seld and all it's consecutive and 
%divided by the ensemble size to get the autocorrelation for the right side
for i = 0 : num_samples-random_check
    stat_acf2(1,i+1) = sum((data(:,random_check) .* data(:,i+random_check))) / ensemble_size;
end
%flip the previous generated result and concatenate it to get the final
%result
stat_acf2 = cat (2, fliplr(stat_acf2(2:end)), stat_acf2);

%polting the two different statistical autocorrelation to compare them so
%that if they are approximatly the same then the process is stationary
%defining a new variable tau to represent the x-axis
tau = (-num_samples+1)*sampling_time :sampling_time: (num_samples-1)*sampling_time;
%ploting the acf of the first column in the left side of the figure
figure ('Name', 'Stationary check');
subplot(1,2,1);
plot(tau, stat_acf);
title ("Statistical ACF for first column");
ylim([-2,20]);
xlabel("Tau(ms)")
ylabel("Autocorrelation");

tau2 = (-num_samples+random_check)*sampling_time :sampling_time: (num_samples-random_check)*sampling_time;
%ploting the acf of the random column in the left right of the figure
subplot(1,2,2);
plot(tau2, stat_acf2);
title ("Statistical ACF for random column");
ylim([-2,20]);
xlabel("Tau(ms)")
ylabel("Autocorrelation");

%Q4: Time Mean and autocorrelation function for one waveform
%1. Time mean:
%choose a random realization to get it's time mean
random_realization = randi([1,ensemble_size]);
%get the mean by getting the sum of the row of the random realization and
%divided by the number of elements which is represented this time by the
%number of samples
time_mean = sum(data(random_realization,:)) / num_samples;
%replicate the mean over the number of samples to be plotted
time_mean = repelem(time_mean,num_samples);

%2. Time Autocorrelation
%firstly create to empty matrix, one to store the shifted data and the
%other to store the result of the time autocorrelation function
shifted = zeros(size(data(random_realization,:)));
time_acf = zeros(size(data(random_realization,:)));

%secondly, according to the chosen line code a random bit will be chosen
%to be put when the waveform is shifted
switch Line_code
    case 1
        %put A or 0 for Unipolar NRZ
        random_bit = randi([0,1]) * A;
    case 2
        %put A or -A for Polar NRZ
        random_bit = (2*(randi([0,1]))-1)*A;
    case 3
        %put A or -A for Polar RZ
        random_bit = (2*(randi([0,1]))-1)*A;
end
%define a new variable to calculate the number of the sample in the bit (from 1 to 7)
flag = 1; 

%loop on the number of samples to multiply the waveform by itself and
%shifted waveform from being shifted by 1 to the number of samples to cover
%it all
for i = 1 : num_samples
    %firstly store the waveform itself to be initially multiply by itself
    if i == 1
        shifted = data(random_realization,:);
        flag = flag+1; %increment the flag after the first sample in the bit is used
    %check if we are using polar return to zero to add the shifted sample
    %firstly equal to zero for three samples then a random value between A
    %or -A for the 4 remaining samples in the bit
    elseif Line_code==3
        if flag == 1 || flag == 2 || flag == 3
            shifted (2:end) = shifted(1: end-1);
            shifted (1) = 0;
        else
            shifted (2:end) = shifted(1: end-1);
            shifted (1) = random_bit;
        end
        flag = flag+1; % increament the flag each sample to determine its position in the bit
    else
        % for unipolar or polar NRZ the shifted bit will be the random bit
        % that has been set before
        shifted (2:end) = shifted(1: end-1);
        shifted (1) = random_bit;
    end
    
    %multiply the waveform by its shifted version then sum and divide it to get the 
    %acf and then store it 
    time_acf(i) = sum((data(random_realization,:) .* shifted(1:end)), 'all') / num_samples;
    
    %check if the bit has been ended and its 7 samples has already shitfted
    %to generate a new random value
    if mod(i,7) == 0
        switch Line_code
            case 1
                random_bit = randi([0,1]) * A;
            case 2
                random_bit = (2*(randi([0,1]))-1)*A;
            case 3
                random_bit = (2*(randi([0,1]))-1)*A;
        end
        flag = 1;
    end
end
%after getting the right sided acf, flip it and concatenate to get the
%final graph
time_acf = cat (2, fliplr(time_acf(2:num_samples)), time_acf);

%Q5: now we can determine if it is ergodic by comparing the figure of the
%time mean with the statistical mean and the time autocorrelation with the
%statistical autocorrelation if both are equal so we can call it an ergodic
%process
figure ('Name', 'Time mean and autocorrelation');
%plotting the time mean in the left half of the figure
subplot(1,2,1);
plot(t, time_mean);
title ("Time Mean");
ylim([-6,6]);
xlabel("Time(ms)")
ylabel("Mean(V)");
%plotting the time autocorrelation funtion in the other half
subplot(1,2,2);
plot(tau, time_acf);
title ("Time Autocorrelation");
ylim([-2,20]);
xlabel("Tau(ms)")
ylabel("Autocorrelation");

%Q6: calculating bandwidth
%to calculate the bandwindth a graph of the PSD has to be plotted to get
%the bandwidth from it
%firstly, to get the x-axis correctly discribe the value of the frequansy
%a variable k has been defined with the size of our acf and then through
%the equation k*fs/(2*num_samples) we can get the right representaion
k = -num_samples + 1: num_samples - 1;
fs = 100;  %Sampling frequency
figure('Name','PSD');
%on the other hand the y-axis is represented by the absolute value of the
%shifted fast fourier transform of the acf to represent the PSD
plot(k*fs/(2*num_samples),abs(fftshift(fft(stat_acf))));
title("PSD");
%limiting the y-axis in unipolar NRZ as it goes to infinity at frequensy = 0
if Line_code == 1
    ylim([0,150]);
end
xlabel("Frequency(Hz)");
ylabel("PSD");