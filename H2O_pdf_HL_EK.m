function [pdf,sigma] = H2O_pdf_HL_EK(x,D0,fitoption,T,Acc,b1,b2,b3)
% x is spatial grid on which you get output; needs regular spacing. 
% DO is surface diffusivity
% T is teperature in K
% Acc is accumulation rate in m/yr
% fitoption gives you option to choose between exponential, lognormal fit - Can be expanded with more options
% b1-b3 are options for the fitting routine - these you can define yourself

dx = x(2)-x(1);

% ----- find firn density, porosity
dz=0.1;
z=0:dz:130;
rho_ice = 0.92; rho_0=0.4;
rho = HL_steady(Acc,T,z,dz,rho_0,rho_ice);  %density using Herron Langway
s = 1-rho/rho_ice;                      %porosity

% ----- calculate ice age
for teller = 1:length(z)
    IceAge(teller) = Trapezoid(rho(1:teller),dz)/(rho_ice*Acc);
end

% ----- find firn close-off
rho_co = 1/( 1/(rho_ice) + T*6.95E-4 - 0.043); % Martinerie et al 1994
z_lid = interp1(rho,z,rho_co-0.012); % find lock-in depth in m
t = floor(interp1(z,IceAge,z_lid)); % find lock-in ice age

% ----- open/closed porosity using Schwander et al.
s_closed = s.*exp(75*(rho/rho_co-1)  ).*(rho<rho_co) + s.*(rho>=rho_co);
s_open = s-s_closed;
s_open_lid = interp1(z,s_open,z_lid);

% ----- find inverse tortuosity:
InvTort = max((s_open-s_open_lid),0);

switch fitoption
    
    case 'exp' % exponential increase;
        % set grain metamorphism parameters
        k=b1;       % a rate parameter
        dt = 0.01/k;    % a characteristic time, divided by 100
        Nt = floor(t/dt);   % number of characteristic times that occur before lock in age
                            % aka number of time steps to travel through
                            % the firn
        
        % find diffusivity as function of time:
        D = interp1(IceAge,InvTort,dt*(1:Nt)); % Inverse tortuosity with time as experienced by a grain
        D = D.*interp1(IceAge,rho,dt*(1:Nt))/rho_ice; % VERY CRUDE correction for strain. this can definitely be improved
        D = D*D0; % scale diffusivity by an input parameter --> makes it very small. Maybe this accounts for units?
        
        % initialize the probability density function
        pdf = zeros(size(x));
        sigma = zeros(1,Nt);
        mu = 0;  % mean x of distribution. Leave this at zero, no reason to change it.
        
        % at each time step add more diffusing water molecules to the open pores
     %  figure;
        for teller =1:Nt
            D_temp = sum(D(teller:end)); % at each time (depth) step, integrate D above
            % should the above be D_temp = sum(D(1:teller)); instead of teller:end? This would
            % integrate the diffusivity down to the depth of the teller,
            % which seems to be what needs to happen here? OR is it from
            % the bottom up because there is more pore space as you move
            % from close off to the surface?
            sigma_temp = sqrt(2*D_temp*dt);  % sigma depends on integrated diffusivity up to that point
            pdf = pdf + dt*k*exp(-k*(teller-0.5)*dt)/sqrt(2*pi)/sigma_temp*exp(-(x-mu).^2/(2*sigma_temp^2));
            sigma(teller) = sigma_temp; % record sigma for each time step
%             hold on
%             plot(x,pdf)
        end
        
        % find center of curve in case you set mu different from zero.
        [dummy, where] = min(abs(x-mu));
                
        % add a delta spike for the molecules that did not diffuse at all:
        pdf(where) = pdf(where) + exp(-k*t)/(x(2)-x(1));
        
%         norm_con = sum(pdf);
%         norm_pdf = pdf/norm_con;
%         figure(2); plot(x,norm_pdf)
%         
%         figure(1);
%         plot(x,pdf)
%         xlabel('Distance on grid')
%         ylabel('Probability of particle at location')
%         title('PDF Exp Distribution')
       
        
%        % for comparison a Gaussian function:
%         D_temp = sum(D);
%         sigma = sqrt(0.5)*sqrt(2*D_temp*dt);
%         pdf_g = 1/sqrt(2*pi)/sigma*exp(-(x-mu).^2/(2*sigma^2));
%         
%         norm_con_g = sum(pdf_g);
%         norm_pdf_g = pdf_g/norm_con_g;
%         
%         figure(2)
%         hold on
%         plot(x,norm_pdf_g)
%         legend('Exponential','Gaussian')
%         
%         figure(1)
%         hold on
%         plot(x,pdf_g)
%         legend('Exponential','Gaussian')
        
    case 'logn' %lognormal distribution
        % set grain metamorphism parameters
        
        dt = 0.01/b1;
        Nt = floor(t/dt);
        
        % find diffusivity as function of time:
        D = interp1(IceAge,InvTort,dt*(1:Nt)); % Inverse tortuosity with time
        D = D.*interp1(IceAge,rho,dt*(1:Nt))/rho_ice; % VERY CRUDE correction for strain, can be improved
        D = D*D0;
        
        % initialize the probability density function
        pdf = zeros(size(x));
        mu = 0;
        
        % define the lognormal distribution we'll be using here
        % note that the width can be set to various values
        % Y = lognpdf(X,MU,SIGMA) returns values at X of the lognormal pdf with 
        % distribution parameters MU and SIGMA. MU and SIGMA are the mean and 
        % standard deviation, respectively, of the associated normal distribution.
        lognpatter = dt*lognpdf(dt*(1:Nt),log(1/b1),b2);

        % at each time step add more diffusing water molecules to the open pores
        for teller =1:Nt
            D_temp = sum(D(teller:end));
            sigma = sqrt(2*D_temp*dt);
            
            pdf = pdf + lognpatter(teller)/sqrt(2*pi)/sigma*exp(-(x-mu).^2/(2*sigma^2));
        end
        
        [dummy, where] = min(abs(x-mu));
        
        % add a delta spike for the molecules that didn't diffuse yet
        pdf(where) = pdf(where) + (1-sum(dx*pdf));
        
%          figure;
%         plot(x,pdf)
%         xlabel('Distance on grid')
%         ylabel('Probability of particle at location')
%         title('PDF Lognormal Distribution')
        
%         % for comparison:
%         D_temp = sum(D);
%         sigma = sqrt(0.5)*sqrt(2*D_temp*dt);
%         pdf_g = 1/sqrt(2*pi)/sigma*exp(-(x-mu).^2/(2*sigma^2));
        
%         hold on
%         plot(x,pdf_g)
%         legend('Lognormal','Gaussian')

    case 'delta' % sigma distribution is a delta function - this is the classical Gaussian diffusion following Johnsen
      
        mu = 0;
        sigma = b1; % first fitting parameter input is the diffusion length
        if sigma>0  % if diffusion length is greater than zero
            pdf = 1/sqrt(2*pi)/sigma*exp(-(x-mu).^2/(2*sigma^2));   % Gaussian filter to convole the isotope profile from Johnsen
                    % represents probability of a molecule ending at a
                    % certain distance away. Larger sigma -> wider Gaussian
        else        % if diffusion length is zero, ie if molecules did not move
            [dummy, where] = min(abs(x-mu));    % find location of middle of grid
            % add a delta spike for the molecules that did not diffuse at all:
            pdf(where) = 1;     % create pdf that is zero everywhere but 1 at middle of grid x 
                                % because zero probability that molecules
                                % moved any distance away from start
        end
        
end
