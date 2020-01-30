disp('never worked - and possibly useless')

mu = 1;

rhs_t = @(t,x) [x(2,:);
    mu*(1-x(1,:).^2).*x(2,:)-x(1,:)];

rhs = @(x) [x(2,:);
    mu*(1-x(1,:).^2).*x(2,:)-x(1,:)];

y_start = [-1.660361182934604   0.685623238242000];
t_span = linspace(0, 6.677,101);

[t_out, y_out] = ode45( rhs_t, t_span, y_start);

omega = 1/(2*pi*t_out(end));

y_in = y_out';
K = [-(length(y_in(1,:))-1)/2:(length(y_in(1,:))-1)/2];
K_shift = ifftshift(K);
K_big = repmat(K_shift,size(y_in,1),1);

y_fft = ifftshift(fft(y_in,[],2));
 
vector_handle_Fourier = @(y) 0*1i*2*pi*K_big*omega.* y - fftshift(fft( rhs(ifft(ifftshift(y),[],2)),[],2));
rhs_Fourier = @(y) 0*1i*2*pi*K_big*omega.* y - [y(2,:);
    mu * (y(2,:)-conv(conv(y(1,:),y(1,:),'same'),y(2,:),'same'))-y(1,:)];





[y_auto,omega_auto] = Newton_Fourier(vf,x, omega);

return



h = 2*pi/21;
t = 0:h:2*pi;
x = [sin(t(1:end-1));cos(t(1:end-1))]*(1+10^-5*rand);
sink = 0*rand(2,1)*10^-4;
vf = @(x) [x(2,:);mu*(1-x(1,:).^2).*x(2,:)-x(1,:)]+sink;

x_approx = x;
for i = 2:length(h)
    x_approx(:,i) = x_approx(:,i-1) + h*vf(x_approx(:,i-1));
end


%y = Newton_timeseries(vf,x_approx,h);


omega = 2*pi+ 10^-3*rand;
[y_auto,omega_auto] = Newton_Fourier(vf,x, omega);
norm(y_auto)

return
%plot(x(1,:),x(2,:),'r*')
%hold on
%plot(y(1,:),y(2,:),'bo')



% in Fourier space, the functions we deal with would be
y = fftshift([fft(x(1,:));fft(x(2,:))]);
y_test = fftshift(fft(x,[],2));
K = [-(length(y(1,:))-1)/2:(length(y(1,:))-1)/2];
% K_minu_1=1./K;
% K_minu_1((length(y(1,:))+1)/2)=0;
rhs_fft = @(y)[2*pi*1i.*K.*y(1,:)-y(2,:)
    2*pi*1i.*K.*y(2,:)+y(1,:)];
reshape_loc = @(y) [y(1:size(x,2)); 
    y(size(x,2)+1:end)];
squeeze_function = @(z) [z(1,:),z(2,:)].';
reshape_rhs = @(y) squeeze_function(rhs_fft(reshape_loc(y)));
Drhs = @(y) [diag(1i.*2*pi*K), -diag(1+0*x(1,:))
     diag(1+0*x(1,:)),diag(1i.*2*pi*K)];
reshape_Drhs = @(y) Drhs(reshape_loc(y));
y = Newton_handle(reshape_rhs,squeeze_function(y),reshape_Drhs);
% converges to the right thing IF RESHAPE_LOC IS WRONG

% in Fourier space, the functions we deal with would be
K_big = repmat(K,size(x,1),1);
rhs = @(y) 2*pi*1i*K_big.*y - fft( vf(ifft(y,[],2)),[],2);
reshape_loc = @(y) [y(1:size(x,2)).'; 
    y(size(x,2)+1:end).'];
squeeze_function = @(z) [z(1,:),z(2,:)].';
reshape_rhs = @(y) squeeze_function(rhs_fft(reshape_loc(y)));

Drhs = @(y) [diag(1i.*2*pi*K), -diag(1+0*x(1,:))
     diag(1+0*x(1,:)),diag(1i.*2*pi*K)];
reshape_Drhs = @(y) Drhs(reshape_loc(y));
y = Newton_handle(reshape_rhs,squeeze_function(y),[]);






