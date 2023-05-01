% Wiener Filter for Audio File by Imported Audio and Generated Noises
% by: M. Alifuddin Akbar [https://github.com/Shemesty10] 

clear all
clc

[Amplitudo_Sampling,Frekuensi_Sampling] = audioread("Recorded.wav");
% Amplitudo_Sampling(:,2) = Amplitudo_Sampling(:,1);

Random_Noise = randn(size(Amplitudo_Sampling))/32;
Audio_Dengan_Noise = Amplitudo_Sampling + Random_Noise;

% sound(Audio_Dengan_Noise,Frekuensi_Sampling)

Audio_Time = Plot_Audio(Amplitudo_Sampling,Frekuensi_Sampling,Random_Noise,Audio_Dengan_Noise);

Plot_Spectogram(Amplitudo_Sampling(:,1),Frekuensi_Sampling)

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
Parameter_Filter_H = inv(Toep_Ryy)*(Rsy');

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

    Estimasi_Audio(k) = Parameter_Filter_H(1)*Audio_Dengan_Noise_PlusNol(k)...
                   + Parameter_Filter_H(2)*Audio_Dengan_Noise_PlusNol(k1)...
                   + Parameter_Filter_H(3)*Audio_Dengan_Noise_PlusNol(k2);
end

clear k
clear k1
clear k2
Audio_Terfilter = Estimasi_Audio(3:length(Estimasi_Audio));

Plot_Spectogram(Audio_Terfilter,Frekuensi_Sampling)
Plot_Spectogram(Audio_Dengan_Noise(:,1),Frekuensi_Sampling)
%% End of Filter Program 

%Function untuk membuat Spektogram
%Untuk Input Amplitudo Sampling Stereo Tolong Pilih Salah Satu Kolom
function Plot_Spectogram(Input_Data_Amplitudo,Input_Data_Frequensi)

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
