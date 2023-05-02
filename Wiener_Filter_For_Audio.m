% Wiener Filter for Audio File by Imported Audio and Generated Noises
% by: M. Alifuddin Akbar [https://github.com/Shemesty10] 

clear all
clc
clf

[Amplitudo_Sampling,Frekuensi_Sampling] = audioread("Recorded.wav");
% Amplitudo_Sampling(:,2) = Amplitudo_Sampling(:,1);

Random_Noise = randn(size(Amplitudo_Sampling))/64;
Audio_Dengan_Noise = Amplitudo_Sampling + Random_Noise;

% sound(Audio_Dengan_Noise,Frekuensi_Sampling)

Audio_Time = Plot_Audio(Amplitudo_Sampling,Frekuensi_Sampling,Random_Noise,Audio_Dengan_Noise);

% Plot_Spectogram(Audio_Terfilter,Frekuensi_Sampling)
[SpectogramDb_Noised,Freq_Vec_Noised,Tim_Vec_Noised] = Plot_Spectogram(Audio_Dengan_Noise(:,1),Frekuensi_Sampling);
% [SpectogramDb,Freq_Vec,Tim_Vec] = Plot_Spectogram(Amplitudo_Sampling(:,1),Frekuensi_Sampling);
% figure
[Spectogram_Noise,~,~] = Plot_Spectogram(Random_Noise(:,1),Frekuensi_Sampling);

Spectogram_Bersih = SpectogramDb_Noised - Spectogram_Noise;

%% Pembuatan Wiener Filter (Bukan Recursive)
%Harus memperhitungkan Cross Correlation antara Dengan dan Tanpa Gangguan
%Serta Memperhitungkan Autokorelasi Keluaran dengan Gangguan

%y Merupakan Sinyal yang telah terkena Noise dan 
%s merupakan Sinyal Asli yang belum terkena Noise

for n = 1:3
    m = n-1;
    Ryy(n) = Perhitungan_Cross_Correlation(m,Audio_Dengan_Noise(:,1),Audio_Dengan_Noise(:,1));
end

for n = 1:3
    m = n-1;
    Rsy(n) = Perhitungan_Cross_Correlation(m,Amplitudo_Sampling(:,1),Audio_Dengan_Noise(:,1));
end

clear m
clear n

%Ryy ditampilkan dengan Matriks Toeplitz
Toep_Ryy = toeplitz(Ryy');
Parameter_Filter_H_NonRecursive = inv(Toep_Ryy)*(Rsy');

%Urutan data ditambahkan 2 Data didepannya agar array mulai dari 1 (Nanti
%hasil estimasi akan dikurangi 2 sekuen setelah selesai
Audio_Dengan_Noise_PlusNol = [zeros(2);
                     Audio_Dengan_Noise];

Jumlah_Data = length(Audio_Dengan_Noise_PlusNol);

%Perhitungan Hasil Estimasi dengan Menggunakan Filter
Estimasi_Audio(1) = 0;
Estimasi_Audio(2) = 0;
for k = 3:Jumlah_Data
    k1 = k-1;
    k2 = k-2;

    Estimasi_Audio(k) = Parameter_Filter_H_NonRecursive(1)*Audio_Dengan_Noise_PlusNol(k)...
                   + Parameter_Filter_H_NonRecursive(2)*Audio_Dengan_Noise_PlusNol(k1)...
                   + Parameter_Filter_H_NonRecursive(3)*Audio_Dengan_Noise_PlusNol(k2);
end

clear k
clear k1
clear k2
Audio_Terfilter_Static_Parameter = Estimasi_Audio(3:length(Estimasi_Audio));

%% End of Filter Non Recursive Program 

%% Pembuatan Recursive Wiener Filter
%Inisiasi Covariance P
Covariance_P(:,:,1) = 0.001*eye(3);
Covariance_P(:,:,2) = 0.001*eye(3);
Forgetting_Factor = 0.85;

% Inisiasi Parameter H
Paremeter_Filter_H_Recursive(:,:,1) = [0 0 0]';
Paremeter_Filter_H_Recursive(:,:,2) = [0 0 0]';

Amplitudo_Sampling_PlusNol = [zeros(2);Amplitudo_Sampling];
Panjang_PlusNol = length(Amplitudo_Sampling_PlusNol);

for n = 3:Panjang_PlusNol
    Audio_Dengan_Noise_PlusNol_Potong(:,:,n) = [Audio_Dengan_Noise_PlusNol(n) Audio_Dengan_Noise_PlusNol(n-1) Audio_Dengan_Noise_PlusNol(n-2)]';

    Covariance_P_num(:,:,n) = Covariance_P(:,:,n-1)*Audio_Dengan_Noise_PlusNol_Potong(:,:,n)*Audio_Dengan_Noise_PlusNol_Potong(:,:,n)'*Covariance_P(:,:,n-1);
    Covairance_P_den(:,:,n) = (1/Forgetting_Factor) + (Audio_Dengan_Noise_PlusNol_Potong(:,:,n)'...
                               *Covariance_P(:,:,n-1)*Audio_Dengan_Noise_PlusNol_Potong(:,:,n));
    
    Covariance_P(:,:,n) = Covariance_P(:,:,n-1) - ...
                          (Covariance_P_num(:,:,n)/Covairance_P_den(:,:,n));

    %denumerator untuk Covariance_P dan Gain_K sama dan akan di reuse
    Gain_K_num(:,:,n) = Covariance_P(:,:,n-1)*Audio_Dengan_Noise_PlusNol_Potong(:,:,n);
    Gain_K_den(:,:,n) = Covairance_P_den(:,:,n);
    Gain_K(:,:,n) = Gain_K_num(:,:,n)/Gain_K_den(:,:,n);
    
    %Konstanta Filter H dihitung dengan memakai Error sementara
    %Sehingga diperlukan mennghitung error sementara
    for a = 1:3
        Estimasi_Sinyal_inArray(a) = Paremeter_Filter_H_Recursive(a,1,n-1)*...
                            Audio_Dengan_Noise_PlusNol_Potong(n);
    end
    
    Estimasi_Sinyal(n) = sum(Estimasi_Sinyal_inArray,'all');
    Error_Sementara(n) = Amplitudo_Sampling_PlusNol(n) - Estimasi_Sinyal(n);
    
    Paremeter_Filter_H_Recursive(:,:,n) = Paremeter_Filter_H_Recursive(:,:,n-1)...
                                          + Gain_K(:,:,n)*Error_Sementara(n);

    Estimasi_Sinyal_Beneran(n) = Paremeter_Filter_H_Recursive(:,:,n)'*Audio_Dengan_Noise_PlusNol_Potong(:,:,n);
    Error_Estimasi(n) = Amplitudo_Sampling_PlusNol(n) - Estimasi_Sinyal_Beneran(n);

end

clear a

%% End of Filter Recursive Program

%% For_Plotting


%Function untuk membuat Spektogram
%Untuk Input Amplitudo Sampling Stereo Tolong Pilih Salah Satu Kolom
function [Spectogram_decibel,Frequency_Vector,Time_Vector] = Plot_Spectogram(Input_Data_Amplitudo,Input_Data_Frequensi)

    window = hamming(512); % Hamming window of length 512 samples
    overlap = round(0.75 * length(window)); % 75% overlap between windows
    nfft = 1024; % Number of FFT points
    
    [Spectogram_Matrix, Frequency_Vector, Time_Vector] ...
        = spectrogram(Input_Data_Amplitudo, window, overlap, ...
        nfft, Input_Data_Frequensi, 'yaxis');
    
    Spectogram_decibel = 10*log10(abs(Spectogram_Matrix));
    
    figure
    imagesc(Time_Vector,Frequency_Vector,Spectogram_decibel)
    axis xy;
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    colormap(jet);
    colorbar;

end

%Perhitungan Autokorelasi Dengan Lag Tertentu

function [Nilai_Korelasi] = Perhitungan_Cross_Correlation(Lag_Sekian,Data1,Data2)

    [Korelasi,lag] = xcorr(Data1,Data2);
    Urutan_Data = 1;
    Nilai_Korelasi = 0;
    for i = min(lag):1:max(lag)
        if i == Lag_Sekian
            Nilai_Korelasi = Korelasi(Urutan_Data);
        end
    Urutan_Data = Urutan_Data + 1;
    end

end

function [Lama_Audio] = Plot_Audio(Amplitudo_Sampling,Frekuensi_Sampling,Random_Noise,Audio_Dengan_Noise)

    %Frekuensi Sampling menunjukkan bahwa terdapat 48k data dalam sedetik
    %Time Sampling = 1/Frekuensi_Sampling

    Total_Data_Sampling = length(Amplitudo_Sampling);
    Time_Sampling = 1/Frekuensi_Sampling;
    Lama_Audio = Total_Data_Sampling/Frekuensi_Sampling;
    Panjang_Sumbu_x = linspace(0,Lama_Audio,Total_Data_Sampling);

    figure
    subplot(3,2,1)
    plot(Panjang_Sumbu_x,Amplitudo_Sampling(:,1));
    xlim([0 Lama_Audio]);
    xlabel('Time(s)')
    ylabel('Amplitudo')
    title('Suara Rekaman Kiri')

    subplot(3,2,2)
    plot(Panjang_Sumbu_x,Amplitudo_Sampling(:,2));
    xlim([0 Lama_Audio]);
    xlabel('Time(s)')
    ylabel('Amplitudo')
    title('Suara Rekaman Kanan')

    subplot(3,2,3)
    plot(Panjang_Sumbu_x,Random_Noise(:,1));
    xlim([0 Lama_Audio]);
    xlabel('Time(s)')
    ylabel('Amplitudo')
    title('Noise Kiri')

    subplot(3,2,4)
    plot(Panjang_Sumbu_x,Random_Noise(:,2));
    xlim([0 Lama_Audio]);
    xlabel('Time(s)')
    ylabel('Amplitudo')
    title('Noise Kanan')

    subplot(3,2,5)
    plot(Panjang_Sumbu_x,Audio_Dengan_Noise(:,1));
    xlim([0 Lama_Audio]);
    xlabel('Time(s)')
    ylabel('Amplitudo')
    title('Suara yang dikenakan Noise kiri')

    subplot(3,2,6)
    plot(Panjang_Sumbu_x,Audio_Dengan_Noise(:,2));
    xlim([0 Lama_Audio]);
    xlabel('Time(s)')
    ylabel('Amplitudo')
    title('Suara yang dikenakan Noise Kanan')

end