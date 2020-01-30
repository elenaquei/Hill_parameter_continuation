function [x_out,omega_out] = Newton_Fourier(vector_field, x_in, omega)
% function x = Newton_Fourier(vector_field, x_in, omega)
%
% transform x into its Fourier sequence correspondence and solve the
% periodic orbit in Fourier space
%
% INPUT
% vector_field      handle, \dot x = f(x), takes input in R^n and return in
%                   the same size
% x                 approximate periodic orbit, x\in R^(nxm) x(:,j)
% approximate solution at time hj
% h                 time stepping used for x
% period
% OUTPUT
% x_new             better approximation



y_in = (fft(x_in,[],2));
K = [-(length(y_in(1,:))-1)/2:(length(y_in(1,:))-1)/2];
K_big = repmat(K,size(x_in,1),1);
vector_handle_Fourier = @(y) 1i*2*pi*K_big*omega.* y - fft( vector_field(ifft(y,[],2)),[],2);
%phase_cond= @(y) sum(sum(1i*K.*y.*conj(y_in)));


reshape_loc = @(y) reshape(y,size(y_in.')).';
squeeze_function = @(z) reshape(z.',1,[]).';
%omega_and_y = @(omega,y) [omega; squeeze_function(y)];
%omega_selec = @(vec) vec(1);
%y_selec = @(vec) reshape_loc(vec(2:end));



%reshape_vector_handle_Fourier = @(x) [phase_cond(y_selec(x));
%    squeeze_function(vector_handle_Fourier(omega_selec(x),y_selec(x)))];

reshape_vector_handle_Fourier = @(x) squeeze_function(vector_handle_Fourier(reshape_loc(x)));

omega_and_y_out = Newton_handle(reshape_vector_handle_Fourier, squeeze_function(y_in));


x_out = fft(y_selec(omega_and_y_out));
%omega_out = omega_selec(omega_and_y_out);


return 
