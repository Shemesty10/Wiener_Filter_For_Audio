% Wiener Filter for Audio File 
% by: M. Alifuddin Akbar [] 

[Amplitudo_Sampling,Frekuensi_Sampling] = audioread("Recorded.wav");

Random_Noise = randn(size(Amplitudo_Sampling))/10;


%Frekuensi Sampling menunjukkan bahwa terdapat 48k data dalam sedetik
%Time Sampling = 1/Frekuensi_Sampling

Total_Data_Sampling = length(Amplitudo_Sampling);
Time_Sampling = 1/Frekuensi_Sampling;
Lama_Audio = Total_Data_Sampling/Frekuensi_Sampling;

figure
Panjang_Sumbu_x = linspace (0,Lama_Audio,Total_Data_Sampling);
plot(Panjang_Sumbu_x,Amplitudo_Sampling)
xlim([0 Lama_Audio]);
xlabel('Time (s)')
ylabel('Amplitude')


%Function untuk membuat Spektogram
%Untuk Input Amplitudo Sampling Tolong Pilih Salah Satu Kolom

Plot_Spectogram(Amplitudo_Sampling(:,1),Frekuensi_Sampling);
Plot_Spectogram(Random_Noise(:,1),Frekuensi_Sampling);

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