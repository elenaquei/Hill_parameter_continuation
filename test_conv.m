t = 0:0.000001:2*pi;
x = sin(t) + 0*10^-6 *rand(size(t));
n = 10^6;
x_fft = (fftshift(fft(x)));

x_conv = fftshift(fft((ifft(ifftshift(x_fft)).^n)));
semilogy(abs(real(x_conv)))
hold on
semilogy(abs(imag(x_conv)))
