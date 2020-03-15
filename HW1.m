clear; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

% Visualize original data
figure(1)
title('Unfiltered signal in spatial domain')
for j=1:20
Un(:,:,:)=reshape(Undata(j,:),n,n,n);
isosurface(X,Y,Z,abs(Un),1)
axis([-20 20 -20 20 -20 20]), grid on
end
xlabel('x')
ylabel('y')
zlabel('z')

% Find the center frequency by averaging all signals
ave = zeros(n,n,n);
for i=1:20 % iterate through 20 measurements and find the total sum of the 20 measurement
Un(:,:,:)= fftn(reshape(Undata(i,:),n,n,n)); % reshape each measurement and swtich to frequency domain by applying fftn
ave = ave + Un; % add up each measurement
end
ave = abs(fftshift(ave))./20;         % rearrange the order and find average
ave = ave./max(ave(:));               % normalize the average 
[M I] = max(ave(:));                  % find maximum value and linear index
[x_ind, y_ind, z_ind] = ind2sub(size(ave),I); % convert linear index into subscripts to find center frequency location(index)

% Just to visualize the approximate center frequency index in [Kx Ky Kz]
figure (2)
isosurface(Kx,Ky,Kz,abs(ave),0.75)
axis([-20 20 -20 20 -20 20]), grid on, drawnow
title('Average signal in frequency domain with isovalue 0.75')
xlabel('Kx')
ylabel('Ky')
zlabel('Kz')

%%  
tau = 0.4; % Filter width
k0 = [Kx(x_ind, y_ind, z_ind) Ky(x_ind, y_ind, z_ind) Kz(x_ind, y_ind, z_ind)]; % Frequency of interest-center frequency
filter = exp(-tau * ((Kx-k0(1)).^2 + (Ky - k0(2)).^2 + (Kz - k0(3)).^2)); % Define the filter
filter = ifftshift(filter);                        % unshift the order of the filter

figure(3)
title('20 Marbles positions in spatial domain with isovalue 0.75 ')
xlabel('x')
ylabel('y')
zlabel('z')
marble_index = zeros(20,3);        % 20 rows for 20 measurement and 3 columns for x y z coordinates
for i = 1:20
   un(:,:,:) = reshape(Undata(i,:),n,n,n);
   unt = fftn(un);
   unft = filter.*unt;                             % Apply filter to the signal in frequency domain
   unt = ifftn(unft);                              % Transfer the signal back to spatial domain
   unt = abs(unt)./max(unt(:)); 
   [Max, Ind] = max(unt(:));
   [max_x, max_y, max_z] = ind2sub(size(unt),Ind); 
   marble_index(i,:) = [X(max_x, max_y, max_z) Y(max_x, max_y, max_z) Z(max_x, max_y, max_z)];
   isosurface(X,Y,Z,abs(unt),0.75)
   axis([-20 20 -20 20 -20 20]), grid on, drawnow
end
figure(4)
plot3(marble_index(:,1), marble_index(:,2), marble_index(:,3),'k-');
axis([-20 20 -20 20 -20 20]), grid on
hold on
plot3(marble_index(20,1), marble_index(20,2), marble_index(20,3), 'ro')
title('Marble trajectory with the 20th marble in spatial domain')
xlabel('x')
ylabel('y')
zlabel('z')