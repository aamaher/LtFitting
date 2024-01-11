function [y_fit,amp,tau,chi_out,err_status] = mydeconv(y, t_sig, H, t_IRF,param)

%%% Defining zero vectors for outputs
    tau = zeros(param.order,1);
    amp = zeros(param.order,1);
    sig = zeros(param.order,1);
    err_status = 0;
%%% Defining some parameters
    fontsz = 9;
    dt = diff(t_IRF(1:2));
    N = length(y); 
    f_max =1/dt;
    df=1/(N*dt);
    f = (-f_max/2:df:f_max/2-df)';
%%% Rescaling
    if (param.laser == "SP")
        H = H-mean(H(t_IRF-t_IRF(1)<3.5));
        DarkNoise = mean(y(t_sig-t_IRF(1)<3.5));
        y_n_std = std(y(t_sig-t_IRF(1)<3.5));
    elseif (param.laser == "NKT")
        H = H-mean(H(t_IRF-t_IRF(1)<3.2));
        DarkNoise = mean(y(t_sig-t_IRF(1)<3.2));
        y_n_std = std(y(t_sig-t_IRF(1)<3.2));
    end 
    
    H = H/sum(H(:));
    y = masking(y,4);     %%% To avoid zeros, putting 1 instead of zero (Lackwiz page ...)
    y(y<1) = 1;
    y_OG = y;    
    if (t_sig(end) ~= t_IRF(end) || t_sig(1) ~= t_IRF(1))
        t_end_temp   =  min(t_sig(end),t_IRF(end));
        t_start_temp =  max(t_sig(1),t_IRF(1));
        t_temp  = (t_start_temp:dt:t_end_temp)';
        t_temp_OG  = (t_start_temp:diff(t_sig(1:2)):t_end_temp)';
        % t_temp  = (t_IRF(1):t_sig(2)-t_sig(1):(t_IRF(end)+t_sig(2)-t_sig(1)))';
        y_OG = interp1(t_sig, y, t_temp_OG, 'linear');
        y = interp1(t_sig, y, t_temp, 'linear');
        H = interp1(t_IRF, H, t_temp, 'linear');
        H = H/sum(H(:));
        t_IRF = t_temp;
        t_sig = t_temp_OG;
    else
        y = interp1(t_sig, y, t_IRF, 'linear');
    end
    

%%%% After this point, both the signal in time axes are all aligned
    t_sig = t_sig-t_sig(1);
    t_IRF = t_IRF-t_IRF(1);
    t = t_IRF;

%%% Calculating Scaling parameters using a mono exponential
    temp_e = myconv(H,exp(-t/2.5),0);     %%% Using tau = 4 ns
    tshift = dt*(find(abs(y - max(y)) == 0) - find(abs(temp_e - max(temp_e)) == 0));
    if (length(tshift)>1)
        tshift = tshift(1);
    end
    y_max = max(y);
    y_min = min(y);
    temp_min = min(temp_e);
    temp_max = max(temp_e);
    A_ratio = (y_max - y_min) / (temp_max - temp_min);
    % Debug_plots();

%%% Setting the cost function
    [problem,IC] = get_problem();   %%% IC here for debugging, by using "Costfn(IC)"
    if (param.global == 1)
        rng default % For reproducibility
        gs = GlobalSearch;
        [B,Cost]   = run(gs,problem);
    else
        [B,Cost] = fmincon(problem);
    end
    y_fit = output_t(B);
    
    %%%% Sample original points again
    common_time = ismember(round(t_IRF,2),round(t_sig,2));
    y_fit = y_fit(common_time);
    if (param.Method == 2)
        chi_vec = ((y_fit - y_OG).^2)./noise_std(y_fit);
    else
        chi_vec = ((y_fit - y_OG).^2)./y_fit;
    end
    chi_r_vec = chi_vec/length(chi_vec);
    chi = sum(abs(chi_vec));
    chi_r = sum(abs(chi_r_vec));
    [~,idx_max] = max(y_OG);
    chi_vec_clipped = chi_vec(idx_max:end);
    chi_r_clipped = sum(chi_vec_clipped)/length(chi_vec_clipped);
    chi_vec_corrected = chi_vec_clipped;
    chi_vec_corrected(chi_vec_corrected> 15) = [];
    chi_r_corrected = sum(chi_vec_corrected)/length(chi_vec_corrected);
    chi_out = chi_r_clipped;
    param.tau_dist = get_dist(B);

%%% Fill output vectors with fitting results
    if (param.cont == 1)
        for k = 1:param.order
            tau(k) = B(3+(k-1)*3);
            sig(k) = B(4+(k-1)*3);
            amp(k) = B(5+(k-1)*3);
        end
    else
        for k = 1:param.order
            tau(k) = B(3+(k-1)*2);
            amp(k) = B(4+(k-1)*2);
        end
         [tau, idx ] = sort(tau);  %%% Sort B for better viewing
         amp = amp(idx);
    end
%%% Print results and plot if required flags are on
    if (param.DEBUGGING == 1)
        disp(B);
        disp(strcat('Cost: ' , num2str(Cost) , '. chi:', num2str(chi),'. chi_r: ', num2str(chi_r),', Clipped: ', num2str(chi_r_clipped),', Corrected: ', num2str(chi_r_corrected)));
        disp(strcat('Lifetime: ',num2str(tau),', amplitude: ',num2str(amp),', Sigma: ',num2str(sig)));
        if(param.PLOT_FLAG == 1)
            plot_I;
            keyboard    %%% Pauses to see the result before moving to the next one
        end
    end

    function C = Costfn(b)
        y_calc = output_t(b);
        y_meas = y;
        switch param.Method
            case 1  %%% LMS
                C = sum((y_calc - y_meas).^2);
            case 2  %%% Gaussian noise, weighted LMS
                C = sum(((y_calc - y_meas ).^2)./noise_std(y_calc));
            case 3  %%% MLE with poisson prior distribution
                C = 2*sum(y_meas.*log(y_meas./y_calc)-(y_meas-y_calc));
            case 4  %%% Gaussian noise, weighted LMS
                C = sum(abs(((y_calc - y_meas).^2)./y_calc));
        end
        % if (param.prior == 1)
            % C = C  + ...;
        % end
    end

    function plot_I (~,~) 
        D = (y_fit - y_OG)./sqrt(y_fit);
        exp_dist = param.tau_dist;
        
        gcf = figure(7);sgtitle(strcat('Lifetime Fitting, \chi^{2}_{r} = ', num2str(chi_r),' ->', num2str(chi_r_clipped),' ->', num2str(chi_r_corrected)),'fontsize',13,'fontweight','bold');
        set(gcf,'position',[100,-10,700,700])
        ax1 = subplot(5,1,[1 2]);hold on;
        yyaxis left
            plot(t_sig,y_OG,'b','linewidth',2);
            plot(t_sig,y_fit,'k-.','linewidth',1.5);
            set(gca, 'YScale', 'log')
        if (param.order > 1)
            B_test = zeros(size(B));  
            B_test(1) = B(1);
            B_test(2) = B(2);
            for k = 1:param.order
                if (param.cont == 1)
                    B_test(3) = B(3+(k-1)*3);
                    B_test(4) = B(4+(k-1)*3);
                    B_test(5) = B(5+(k-1)*3);
                else
                    B_test(3) = B(3+(k-1)*2);
                    B_test(4) = B(4+(k-1)*2);
                end
                y_temp = output_t(B_test);
                y_temp = y_temp(common_time);
                plot(t_sig,y_temp,'m--','linewidth',1.5);
            end
        end
        grid on;
        ylabel('Intensity','fontweight', 'bold','fontsize', fontsz);            
        
        yyaxis right
            plot(t_IRF,H,'r','linewidth',1.5);

        legend('Measured','Fitted');
        ax2 = subplot(5,1,3);hold on;
        plot(t_sig,D,'k','linewidth',2)
        ylabel('Deviation','fontweight', 'bold','fontsize', fontsz);            
        grid on;
        % ylim([-10,10]);

        ax4 = subplot(5,1,4);hold on; 
        plot(t_sig,chi_vec,'k','linewidth',2)
        ylabel({"\chi^{2}";"Contribution"},'fontweight', 'bold','fontsize', fontsz);            
        grid on;
        xlabel('Time (ns)','fontweight', 'bold','fontsize', fontsz);
        linkaxes([ax1,ax2,ax4],'x')
        % xlim([t(idx_max+5),max(t)])
        xlim([0,max(t)])
        ax5 = subplot(5,1,5);hold on; 
        if (param.cont == 0)
            exp_dist (exp_dist == 0 ) = NaN;
            stem(param.tau_vec,exp_dist,'k','linewidth',2)
        else
            plot(param.tau_vec,param.tau_dist,'k','linewidth',2);
        end
        grid on;
        xlim([0,param.tau_vec(end)]);
        if (param.cont == 0)
            Str_type = 'Discrete';
        elseif (param.cont_type == 1)
            Str_type = 'Gaussian';
        elseif (param.cont_type == 2)
            Str_type = 'Gamma';
        elseif (param.cont_type == 3)
            Str_type = 'Inverse Gamma';
        elseif (param.cont_type == 4)
            Str_type = 'Log-normal';
        end
        ylabel({"Amplitude";Str_type},'fontweight', 'bold','fontsize', fontsz);            
        hx = xlabel('Lifetime (ns)','fontweight', 'bold','fontsize', fontsz);
        

        if(param.PLOT_Tau)   %%% Plot frequency info
            figure(72);hold on
            if (param.cont == 0)
                exp_dist (exp_dist == 0 ) = NaN;
                stem(param.tau_vec,exp_dist,'linewidth',2)
            else
                plot(param.tau_vec,param.tau_dist,'linewidth',2);
            end
            grid on;
            xlim([0,param.tau_vec(end)]);
            ylabel('Amplitude','fontweight', 'bold','fontsize', fontsz);            
            hx = xlabel('Lifetime (ns)','fontweight', 'bold','fontsize', fontsz);
        end
        if(param.PLOT_Freq)   %%% Plot frequency info
            y_W = FT(y,dt);
            H_W = FT(H,dt);
            output_w = @(b) FT(output_t(b),dt);
            figure(8);hold on;
            plot(f,abs(y_W)/max(abs(y_W)),'b','linewidth',2);
            plot(f,abs(H_W)/max(abs(H_W)),'r-.','linewidth',2);
            plot(f,abs(output_w(B))/max(abs(output_w(B))),'k-.','linewidth',2);
            grid on;
            hx = xlabel('Freq (GHz)');
            hy = ylabel('Intensity');
            leg = legend('Measured','IRF','Fitted');
            figure_set_maher
        end
    end

    function exp_dist = get_dist(B)
        tau_vector = param.tau_vec;
        exp_dist = zeros(size(tau_vector));
        if (param.cont == 0)
            for j = 1:param.order
                tau_temp = B(3+2*(j-1));
                idx = find(abs(tau_vector - tau_temp) == min(abs(tau_vector - tau_temp)));
                exp_dist(idx) = exp_dist(idx) + B(4+(j-1)*2);
            end
        else
            for j = 1:param.order
                mu = B(3+(j-1)*3);
                sigma = B(4+(j-1)*3);
                if (sigma == 0)
                    continue;
                end
                if(param.cont_type == 1)    %%%% Gaussian distribution
                    pdf_temp = normpdf(tau_vector,mu,sigma)*(tau_vector(end)-tau_vector(1))/length(tau_vector);
                elseif(param.cont_type == 2)    %%%% Gama distribution
                    pdf_temp = gampdf(tau_vector,mu^2/sigma,sigma/mu)*(tau_vector(end)-tau_vector(1))/length(tau_vector);
                elseif(param.cont_type == 3)    %%%% Inverse Gama distribution
                    pdf_temp = gampdf(1./tau_vector, mu^2 / sigma, sigma/(mu^3)) * (tau_vector(end) - tau_vector(1)) / length(tau_vector);    
                elseif(param.cont_type == 4)    %%%% Log-normal distribution
                    gamma = asinh(sigma/(2*mu));
                    pdf_temp = exp(-((log(1./tau_vector)-log(1/mu))/gamma).^2);    
                end
                exp_dist = exp_dist + B(5+(j-1)*3)*pdf_temp;
            end
        end
    end
    
    function [problem,IC] = get_problem(~,~)
        lb = [-param.t_shift_max, 1 - param.amp_ratio_max];
        ub = [ param.t_shift_max, 1 + param.amp_ratio_max];
        IC = [0,1];        
        Beq = 1;
        Aeq = [0,0];
        if (param.cont == 0)
        %%% Multi-exponential
        %%% Parameters are: [Shift  Amp   tau-1  Amp-1 ... ]
        %%% Parameters are: [b(1)   b(2)  b(3)   b(4)  ... ]
            for k = 1:param.order
                Aeq =[Aeq,0,1];
                lb = [lb,param.lb,0];
                ub = [ub,param.ub,1];
                IC = [IC,k*param.IC,1/param.order];
            end
        elseif (param.cont == 1)
            for k = 1:param.order
                Aeq =[Aeq,0,0,1];
                lb = [lb,param.lb  ,param.sigma_max/1000,0];
                ub = [ub,param.ub  ,param.sigma_max     ,1];
                IC = [IC,k*param.IC,param.sigma_max/2   ,1/param.order];
            end
        end
        Aeq =[Aeq,0];
        lb = [lb,-DarkNoise];
        ub = [ub,param.y_shift_max];
        IC = [IC,0];
        opts = optimoptions('fmincon', 'Algorithm','interior-point','MaxFunEvals', 50000, 'MaxIter', 10000, 'StepTolerance', 1e-10,'Display', 'off');
        problem = createOptimProblem('fmincon','x0',IC,'objective',@Costfn,'lb',lb,'ub',ub,'Aeq', Aeq, 'beq', Beq,'options', opts);
    end
    
    function output = output_t(B)
        %%% Parameters are: [shift Amp  mu   sigma amp  ... ]
        %%% Parameters are: [b(1)  b(2) b(3) b(4)  b(5) ...]
        exp_vector = zeros(size(y,1),1);
        if (param.cont == 1)
            exp_dist = get_dist(B);
            exp_vector = sum(exp_dist .* exp(-t ./ param.tau_vec),2);
        elseif (param.cont == 0)
            for k = 1:param.order
                exp_vector = exp_vector + B(4+(k-1)*2)*exp(-t/B(3+(k-1)*2));
            end
        end
        exp_vector(isnan(exp_vector))=0;    %%% To avoid any errors when first point is NaN at 0 lifetime
        output  = B(2)*A_ratio*myconv(H,exp_vector,tshift+B(1))+DarkNoise+B(end);
    end 
    
    function Iout = masking (Iin,neighborhood_size )
        threshold = 1;
        mask = Iin <= threshold;
        Iout = Iin;
        for i = 1:numel(Iin)
            if mask(i)
                start_idx = max(1, i - floor(neighborhood_size / 2));
                end_idx = min(numel(Iin), i + floor(neighborhood_size / 2));
                average_value = mean(Iin(start_idx:end_idx));
                Iout(i) = average_value;
            end
        end
    end
    
    function N = noise_std(I)
        RON = 7;
        Dark_std = DarkNoise;
        N = 4*abs(I) + RON + Dark_std; %%% Gaussian noise
    end

    function C = myconv(A, B,t_shift)   % perform 1D Linear discrete convolution followed by clipping
        N = length(A);
        %%% Zero Padding
        N_ZP = 1024;
        A_ZP = cat(1,A,zeros(N_ZP-N,1));
        B_ZP = cat(1,B,zeros(N_ZP-N,1));
        f_max =1/dt;
        df_ZP=1/(N_ZP*dt);
        f_ZP = (-f_max/2:df_ZP:f_max/2-df_ZP)';
        %%% FFT --> shift and Conv --> IFFT
        A_W = FT(A_ZP,dt);
        B_W = FT(B_ZP,dt);
        try
            C_W = A_W.*B_W.*exp(-1j*2*pi*f_ZP*t_shift);
        catch
            if (param.DEBUGGING == 1)
                keyboard;
            else
                err_status = 1;
                return;
            end
        end
        C = abs(IFT(C_W,1/(N*dt)));
        %%% Clipping the zero padded signal
        C = C(1:N);
    end
    
    function G = FT(g, dt)  % 1D DFT where delta is sample width
        G = fftshift(fft(fftshift(g))) * dt;
    end
    
    function g = IFT(G, df) % 1D DIFT where df is sample with in frequency domain (1/N/dt)
        g = (ifft(ifftshift(G)) * length(G) * df);
    
    end

    function Debug_plots(~,~) 
        figure(1);hold on;
        plot(t,A_ratio*myconv(H,exp(-t/2),-1  ),'b');
        plot(t,A_ratio*myconv(H,exp(-t/2),-0.9),'b');
        plot(t,A_ratio*myconv(H,exp(-t/2),-0.8),'b');
        plot(t,A_ratio*myconv(H,exp(-t/2),-0.7),'b');
        plot(t,A_ratio*myconv(H,exp(-t/2),-0.6),'b');
        plot(t,A_ratio*myconv(H,exp(-t/2),-0.5),'b');
        plot(t,A_ratio*myconv(H,exp(-t/2),-0.4),'b');
        plot(t,A_ratio*myconv(H,exp(-t/2),-0.3),'b');
        plot(t,A_ratio*myconv(H,exp(-t/2),-0.2),'b');
        plot(t,A_ratio*myconv(H,exp(-t/2),-0.1),'b');
        plot(t,A_ratio*myconv(H,exp(-t/2), 0  ),'k','linewidth',2);
        plot(t,A_ratio*myconv(H,exp(-t/2), 0.1),'r');
        plot(t,A_ratio*myconv(H,exp(-t/2), 0.2),'r');
        plot(t,A_ratio*myconv(H,exp(-t/2), 0.3),'r');
        plot(t,A_ratio*myconv(H,exp(-t/2), 0.4),'r');
        plot(t,A_ratio*myconv(H,exp(-t/2), 0.5),'r');
        plot(t,A_ratio*myconv(H,exp(-t/2), 0.6),'r'); 
        plot(t,A_ratio*myconv(H,exp(-t/2), 1  ),'r');
        keyboard 
    end

end 
