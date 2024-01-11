if (param.SHOW_ONLY == 1)
    if (isfile(filename_saving_results))
        load(filename_saving_results,'I_avg','Peak_Img','param','err_seg','chi_seg','tau_map','err_map','amp_map','chi_map','HorizontalPixels','VerticalPixels','I_fit_seg','I_seg_actual','tau_seg','tau_dist_seg','time','Obj_Mask','-mat')
    else
        disp('Cannot find .mat file!')
        status = 5;
        return;
    end
end

%%% Main file name saving template
fig_name_template = strcat('RESULTS_EXP',num2str(param.order),str_cont,'_',filename_saving,'_M',num2str(param.Method));

%%% Check whether in folder or not
if (param_plot.SaveInFolder == 0)
    fig_name_template = strcat ('Results/',fig_name_template);
elseif (param_plot.SaveInFolder == 1)
    folder1 = strcat('Results/',filename_saving_temp(end-1));
    folder2 = strcat('Results/',filename_saving_temp(end-1),'/',filename_saving_temp(end));
    fig_name_template = strcat (folder2,'/',fig_name_template);    
    if (~isfolder(folder1))
        mkdir(char(folder1));
    end
    if (~isfolder(folder2))
        mkdir(char(folder2));
    end
end

%%% Initialize Variables
tau_plot = tau_map;
au_plot(tau_plot == 0) = NaN;
if (param.order == 1)
    [V,E] = histcounts(tau_plot,param_plot.hbins_tau);
    W = V/sum(V);
    dE = diff(E(1:2));
    tau_mean = sum((E(1:end-1)+dE/2).*W);
    tau_std = sqrt(sum(W.*(E(1:end-1) + dE/2 - tau_mean).^2));
else
    E = param_plot.hbins_tau;
    [V, ~] = histwv(tau_plot, amp_map, 0, E(end), length(param_plot.hbins_tau));
    W  = V/sum(V);
    dE = diff(E(1:2));
    tau_mean = sum((E+dE/2).*W);
    tau_std = sqrt(sum(W.*(E + dE/2 - tau_mean).^2));
end

hbins_chi = 0:0.1:100;
chi_avg = mean(chi_map(chi_map>0));
chi_map (chi_map == 0) = NaN;
[V_chi,E_chi] = histcounts(chi_map,hbins_chi);
dE_chi = diff(E_chi(1:2));
W_chi = V_chi/sum(V_chi);
chi_mean = sum(E_chi(1:end-1).*W_chi);
chi_std = sqrt(sum(W_chi.*(E_chi(1:end-1) + dE_chi/2 - chi_mean).^2));
  
if (param_plot.plot_mix == 1)
    gcf = figure(1);
    clf(gcf);    
    if (param.order == 1)
        ax1 = axes(gcf,'InnerPosition',[0.1300 0.1100 0.72 0.8150] );
        subplot(1,3,2);hold on;
        histogram(tau_plot,param_plot.hbins_tau);
        xlabel('Lifetime (ns)','fontsize', 15,'fontweight', 'bold');
        ylabel('Count','fontsize', 15,'fontweight', 'bold');
        title(['Mean = ' num2str(round(tau_mean,2)) ', Std = ' num2str(round(tau_std,3))],'fontsize', 14,'fontweight', 'bold')
        xline(tau_mean,'Color','red','LineWidth',1.5)
        xlim([max(0,tau_mean-10*tau_std),tau_mean+10*tau_std]);
        pbaspect([1 1 1]);
        IC = [max(V),tau_mean,tau_std];
        if (param_plot.Fit_Histogram == 1)
            pdf_temp = normpdf(E,IC(2),IC(3));
            pdf_temp = IC(1)*pdf_temp/max(pdf_temp);
            plot(E+dE/2,pdf_temp,'r--','LineWidth',2)
        end
        ax2 = subplot(1,3,1);
        imagesc(ax2,HorizontalPixels, VerticalPixels, tau_plot);
        colormap(ax2,'turbo')
        colorbar;
        clim([max(0,tau_mean-3*tau_std),tau_mean+3*tau_std]);
        title('Lifetime (ns)','fontsize', 14,'fontweight', 'bold');
        pbaspect([1 1 1]);
        ax3 = subplot(1,3,3);
        imagesc(HorizontalPixels, VerticalPixels, chi_map);
        colormap(ax3,'turbo')
        title(['\chi^2_{r} (Mean: ' num2str(round(chi_mean,1)) ', Std:' num2str(round(chi_std,1)) ')'],'fontsize', 14,'fontweight', 'bold');
        clim([max(0,chi_mean-1*chi_std),chi_mean+1*chi_std])
        colorbar;
        pbaspect([1 1 1]);
        set(gcf,'position',[10,100,1200,400])

    elseif (param.order > 1)
        colormap jet
        IC = [];
        for i = 1:param.order
            subplot(param.order,3,2+3*(i-1));hold on;
                % h = histogram([tau_plot(:,:,i)],param_plot.hbins_tau);
                [V_temp, ~] = histwv(tau_plot(:,:,i), amp_map(:,:,i), 0, 10, length(param_plot.hbins_tau));
                bar(param_plot.hbins_tau,V_temp,'k','linewidth',2)
                E_temp = param_plot.hbins_tau;
                W_temp = V_temp/sum(V_temp);
                tau_mean_temp = sum(E_temp.*W_temp);
                tau_std_temp  = sqrt(sum(W_temp.*(E_temp + dE/2 - tau_mean_temp).^2));
                xlabel('Lifetime (ns)','fontsize', 13,'fontweight', 'bold');
                ylabel('Count');
                title(['Mean = ' num2str(round(tau_mean_temp,1)) ',Std: ' num2str(round(tau_std_temp,1))])
                xline(tau_mean_temp,'Color','red','LineWidth',1.5)
                xlim([max(0,tau_mean_temp-3*tau_std_temp),tau_mean_temp+3*tau_std_temp]);
                pbaspect([1 1 1]);
                IC = [IC max(V_temp),tau_mean_temp,tau_std_temp];
                subplot(param.order,3,1+3*(i-1));
                imagesc(HorizontalPixels, VerticalPixels, tau_plot(:,:,i));
                colorbar;
                % xlabel('Horizontal Pixels','fontsize', 15,'fontweight', 'bold');
                % ylabel('Vertical Pixels','fontsize', 15,'fontweight', 'bold');
                clim([max(0,tau_mean_temp-3*tau_std_temp),tau_mean_temp+3*tau_std_temp]);
                title(strcat('lifetime ',num2str(i),' (ns)'));
                set(gcf,'position',[100,100,1200,400*param.order])
                pbaspect([1 1 1]);
        end                      
        subplot(param.order,3,3);hold on    
        %%% Weighted histogram here    
        
        [V, ~] = histwv(tau_plot, amp_map, 0, 10, length(param_plot.hbins_tau));
        bar(param_plot.hbins_tau,V,'k','linewidth',2)

        if (param_plot.Fit_Histogram == 1)
            tau_vector = param_plot.hbins_tau;
            tau_vector(end) = [];
            fit_hist = FitHistogram(V,param_plot,param,IC);
            plot(tau_vector+dE/2,fit_hist,'r--','linewidth',2);   
        end
        xlabel('Lifetime (ns)','fontsize', 13,'fontweight', 'bold');
        ylabel('Count');
        title(['Mean = ' num2str(round(tau_mean,1)) ',Std: ' num2str(round(tau_std,1))])
        xline(tau_mean,'Color','red','LineWidth',1.5)
        xlim([max(0,tau_mean-3*tau_std),tau_mean+3*tau_std]);
        pbaspect([1 1 1]);
        subplot(param.order,3,6)
        imagesc(HorizontalPixels, VerticalPixels, chi_map);
        colorbar;
        % hx = xlabel('Horizontal Pixels');
        % hy = ylabel('Vertical Pixels');    
        clim([max(0,chi_mean-1*chi_std),chi_mean+1*chi_std])
        title(['\chi^2_{r} (Mean: ' num2str(round(chi_mean,1)) ', Std:' num2str(round(chi_std,1)) ')'],'fontsize', 15,'fontweight', 'bold');
        % sgtitle('Lifetime with Bi-exponential Fitting' ,'Fontsize',20);
        pbaspect([1 1 1]);
        set(gcf,'position',[100,0,1200,800])
    end
    fig_name = get_file_name(fig_name_template,param_plot);
    saveas(figure(1),fig_name)
    close(figure(1));
end

if (param_plot.plot_chi == 1)
    gcf = figure(1);
    clf(gcf);
    ax = axes(gcf);
    imagesc(ax,HorizontalPixels, VerticalPixels, chi_map);
    colormap jet
    colorbar;
    clim(ax,[max(0,chi_mean-1*chi_std),chi_mean+1*chi_std])
    title(ax,['\chi^2_{r} (Mean: ' num2str(round(chi_mean,1)) ', Std:' num2str(round(chi_std,1)) ')'],'fontsize', 15,'fontweight', 'bold');
    fig_name = get_file_name(fig_name_template,param_plot);
    set(ax,'XColor', 'none','YColor','none')
    
    h2 = colorbar(ax);
    h2.Title.String = '\chi^2';
    h2.Location = 'eastoutside';
    h2.FontSize = 14;
    saveas(figure(1),fig_name)
    close(figure(1));
end

if (param_plot.plot_hist == 1)      %%% Histogram plot only
    gcf = figure(1);hold on;
    if (param_plot.Stack_figure ~= 1)
        clf(gcf);
    end
    E_temp = param_plot.hbins_tau;
    if (param.order > 1)
        [V_temp, ~] = histwv(tau_plot, amp_map, 0, E_temp(end), length(param_plot.hbins_tau));
        W_temp = V_temp/sum(V_temp);
        tau_mean = sum(W_temp.*E_temp);
        tau_std  = sqrt(sum(W_temp.*(E_temp + dE/2 - tau_mean).^2));
        bar(E_temp,V_temp,'linewidth',2);
    else
        E_temp(1) = [];
        bar(E_temp,V,'linewidth',2);
    end
    xlabel('Lifetime (ns)','fontsize', 15,'fontweight', 'bold');
    ylabel('Count','fontsize', 15,'fontweight', 'bold');
    title(['Mean = ' num2str(round(tau_mean,3)) ', Std = ' num2str(round(tau_std,3))],'fontsize', 14,'fontweight', 'bold')
    xline(tau_mean,'Color','red','LineWidth',1.5)
    xlim([max(0,tau_mean-0.5),tau_mean+0.5]);
    if (param_plot.Fit_Histogram == 1)
        IC = [];
        if (param.order > 1)
            for i = 1:param.order
                [V_temp,E_temp] = histcounts(tau_plot(:,:,i),param_plot.hbins_tau);
                V_temp(1) = 0;
                W_temp = V_temp/sum(V_temp);
                tau_mean_temp = sum(E_temp(1:end-1).*W_temp);
                tau_std_temp  = sqrt(sum(W_temp.*(E_temp(1:end-1) + dE/2 - tau_mean_temp).^2));
                IC = [IC max(V_temp),tau_mean_temp,tau_std_temp];
            end
            tau_vector = param_plot.hbins_tau;
            fit_hist = FitHistogram(V,param_plot,param,IC);            
            plot(tau_vector,fit_hist,'r--','linewidth',2);   
        elseif (param.order == 1)
            IC = [max(V),tau_mean,tau_std];
            fit_hist = normpdf(E(2:end)+dE/2,tau_mean,tau_std);
            fit_hist = max(V)*fit_hist./max(fit_hist(:));
            plot(E(2:end)+dE,fit_hist,'r--','LineWidth',1.5)
        end
    end
    if (param_plot.Stack_figure ~= 1)
        fig_name = get_file_name(fig_name_template,param_plot);
        saveas(figure(1),fig_name)
        close(figure(1));
    end
end

if (param_plot.plot_PeakInt == 1)   %%% Peak Intensity
    gcf = figure(1);hold on;
    gcf.Position = [50 50 980 800];
    ax1_1 = axes(gcf);
    % ax1_1 = axes(gcf,'InnerPosition',[0.1300 0.1100 0.72 0.8150] );
    ax1_2 = copyobj(ax1_1,gcf);
    imagesc(ax1_1,HorizontalPixels, VerticalPixels, Peak_Img);
    colormap(ax1_1,gray)
    Scalebar_length = param_plot.scale_bar;
    PixelSize = 11.22/param_plot.Magnification;
    Size_pixels = round(Scalebar_length/PixelSize);    
      
    if (param_plot.Crop_plot == 2)
        Mask = (tau_plot > 0);
        xmin = find(sum(Mask,1), 1, 'first');
        ymin = find(sum(Mask,2), 1, 'first');
        xmax = find(sum(Mask,1), 1, 'last');
        ymax = find(sum(Mask,2), 1, 'last');

        xmin = max(0,xmin-25);
        xmax = min(size(Mask,2),xmax+25);
        ymin = max(0,ymin-25);
        ymax = min(size(Mask,1),ymax+25);

        % x_vector = [xmax-Size_pixels-40 xmax-40];   % for plotting
        x_vector = [xmin+ 20 xmin+Size_pixels+20];   % for plotting
        y_vector = [ymax-20 ymax-20];   % for plotting
        line(ax1_2,x_vector,y_vector,'LineWidth',param_plot.scale_bar_Linewidth,'Color',[0.99 0.99 0.99]);  % plot line
        xlim(ax1_1,[xmin,xmax]);
        ylim(ax1_1,[ymin,ymax]);
        pbaspect(ax1_1,[xmax-xmin ymax-ymin 1])
        pbaspect(ax1_2,[xmax-xmin ymax-ymin 1])        
    elseif (param_plot.Crop_plot == 0)
        xmax = 960;
        ymax = 608;
        xmin = 0;
        x_vector = [xmin+20 xmin+Size_pixels+20];   % for plotting
        y_vector = [ymax-20 ymax-20];   % for plotting
        line(ax1_2,x_vector,y_vector,'LineWidth',param_plot.scale_bar_Linewidth,'Color',[0.99 0.99 0.99]);  % plot line
    end
    ax1_2.UserData   = linkprop([ax1_1,ax1_2],{'x','y','Position','InnerPosition','xtick','ytick','ydir','xdir','xlim','ylim'});
    set(gcf.Children,'XColor', 'none','YColor','none')
    set(gcf.Children,'Visible', 'off')

%%% Saving
    fig_name = get_file_name(fig_name_template,param_plot);
    saveas(figure(1),fig_name)
    close(figure(1));
end

if (param_plot.plot_PeakLT == 1)    %%% Peak Intensity + Lifetime
    
    if (param.order >1 )
        tau_plot = tau_map(:,:,1);
    end
    gcf = figure(1);hold on;
    % ax2_1 = axes(gcf,'InnerPosition',[0.1300 0.1100 0.72 0.8150] );
    ax2_1 = axes(gcf);
    gcf.Position = [800 50 980 800];
    ax2_2 = copyobj(ax2_1,gcf);
    ax2_3 = copyobj(ax2_1,gcf);
    Image = Peak_Img/max(Peak_Img(:));
    % Image = imsharpen(Image,'Radius',2,'Amount',2);
    % Image = Image/max(Image(:));
    % Image = imbilatfilt(Image,0.5,1.5);
    imagesc(ax2_1,Image/max(Image(:))); %% Plot the image
    clim(ax2_1,[0.02,0.8]);
    imagesc(ax2_2,tau_plot, 'AlphaData', (tau_plot>0));
    % clim(ax2_2,[max(0,tau_mean-3*tau_std),tau_mean+3*tau_std]);     %%% Clim of lifetime
    clim(ax2_2,[1.5,3.5]);     %%% Clim of lifetime
    colormap(ax2_1,'gray')
    colormap(ax2_2,jet(30))
    Scalebar_length = param_plot.scale_bar;
    PixelSize = 11.22/param_plot.Magnification;
    Size_pixels = Scalebar_length/PixelSize;    

    if (param_plot.Crop_plot == 2)
        Mask = (tau_plot(:,:,1) > 0);
        xmin = find(sum(Mask,1), 1, 'first');
        ymin = find(sum(Mask,2), 1, 'first');
        xmax = find(sum(Mask,1), 1, 'last');
        ymax = find(sum(Mask,2), 1, 'last');
        xmin = max(0,xmin-50);
        xmax = min(size(Mask,2),xmax+50);
        ymin = max(0,ymin-50);
        ymax = min(size(Mask,1),ymax+50);
        % x_vector = [xmax-Size_pixels-40 xmax-40];   % for plotting
        x_vector = [xmin+ 20 xmin+Size_pixels+20];   % for plotting
        y_vector = [ymax-20 ymax-20];   % for plotting
        line(ax2_3,x_vector,y_vector,'LineWidth',param_plot.scale_bar_Linewidth,'Color',[0.99 0.99 0.99]);  % plot line
        xlim(ax2_1,[xmin,xmax]);
        ylim(ax2_1,[ymin,ymax]);
        pbaspect(ax2_1,[xmax-xmin ymax-ymin 1])
        pbaspect(ax2_3,[xmax-xmin ymax-ymin 1])        
        pbaspect(ax2_2,[xmax-xmin ymax-ymin 1])        
   
    elseif (param_plot.Crop_plot == 0)
        xmax = 960;
        ymax = 608;  
        xmin = 0;
        % x_vector = [xmax-Size_pixels-20 xmax-20];   % for plotting
        x_vector = [xmin+ 20 xmin+Size_pixels+20];   % for plotting
        y_vector = [ymax-20 ymax-20];   % for plotting
        line(ax2_3,x_vector,y_vector,'LineWidth',param_plot.scale_bar_Linewidth,'Color',[0.99 0.99 0.99]);  % plot line
    end
    ax2_2.UserData = linkprop([ax2_1,ax2_2],{'x','y','Position','InnerPosition','xtick','ytick','ydir','xdir','xlim','ylim'});
    ax2_3.UserData = linkprop([ax2_1,ax2_3],{'x','y','Position','InnerPosition','xtick','ytick','ydir','xdir','xlim','ylim'});
    set(gcf.Children,'XColor', 'none','YColor','none')
    set(gcf.Children,'Visible', 'off')

%%% Colorbar settings
    h2 = colorbar(ax2_2);
    h2.Title.String = '\tau (ns)';
    h2.Location = 'eastoutside';
    h2.FontSize = param_plot.fontsz;
    % h2.Position(1) = h1.Position(1) + 0.1;
    % h2.AxisLocation = 'out';
%%% Saving
    fig_name = get_file_name(fig_name_template,param_plot);
    saveas(figure(1),fig_name)
    close(figure(1));
end

function Manual_Plot_Crop(ax1_1,ax1_2,ax1_3,ax2_1,ax2_2,ax2_3,ax2_4,L1,L2,T1,T2,Size_pixels,Peak_Img,h1,h2,param_plot)
    %%%% Script for using it
    % elseif (param_plot.Crop_plot == 1)   %%% Manual Plot Crop
    %     xmax = 960;
    %     ymax = 608;  
    %     xmin = 0;
    %     x_vector = [xmin+ 20 xmin+Size_pixels+20];   % for plotting
    %     y_vector = [ymax-20 ymax-20];   % for plotting
    %     L1 = line(ax1_2,x_vector,y_vector,'LineWidth',param_plot.scale_bar_Linewidth,'Color',[0.99 0.99 0.99]);  % plot line
    %     L2 = line(ax2_3,x_vector,y_vector,'LineWidth',param_plot.scale_bar_Linewidth,'Color',[0.99 0.99 0.99]);  % plot line
    %     T1 = text(ax2_4,xmax-Size_pixels,ymax-50,bartxt,'color',[0.99 0.99 0.99],'FontSize',14);
    %     T2 = text(ax1_3,xmax-Size_pixels,ymax-50,bartxt,'color',[0.99 0.99 0.99],'FontSize',14);
    %     Manual_Plot_Crop(ax1_1,ax1_2,ax1_3,ax2_1,ax2_2,ax2_3,ax2_4,L1,L2,T1,T2,Size_pixels,Peak_Img,h1,h2,param_plot)
    fig = uifigure;
    fig.Position = [50, 50, 350, 400];

    X = ax1_1.XLim(1);
    Y = ax1_1.YLim(1);
    L = ax1_1.XLim(2) - X;
    W = ax1_1.YLim(2) - Y;
    CB1_P1 = h1.Position(1);
    CB2_P1 = h1.Position(2);
    h = 0;
    w = 120;
    xmin = X;
    xmax = X+L;
    ymin = Y;
    ymax = Y+W;

    % Create a slider
    slider_L = uicontrol('Parent',fig,'Style', 'slider', 'Min', 0, 'Max', size(Peak_Img,2), 'Value', L, ...
        'Position', [w, h+200, 100, 20], 'Callback', @updateImage);
    slider_W = uicontrol('Parent',fig,'Style', 'slider', 'Min', 0, 'Max', size(Peak_Img,1), 'Value', W, ...
        'Position', [w, h+230, 100, 20], 'Callback', @updateImage);
    slider_Y = uicontrol('Parent',fig,'Style', 'slider', 'Min', 0, 'Max', size(Peak_Img,1), 'Value', Y, ...
        'Position', [w, h+260, 100, 20], 'Callback', @updateImage);
    slider_X = uicontrol('Parent',fig,'Style', 'slider', 'Min', 0, 'Max', size(Peak_Img,2), 'Value', X, ...
        'Position', [w, h+290, 100, 20], 'Callback', @updateImage);
    slider_ftsz = uicontrol('Parent',fig,'Style', 'slider', 'Min', 8, 'Max', 50, 'Value', h1.FontSize, ...
        'Position', [w, h+170, 100, 20], 'Callback', @updateImage);

    txt_L = uicontrol('Parent',fig,'Style', 'text', 'Position', [w-100, h+200, 50, 20], ...
        'String', 'Length');
    txt_W = uicontrol('Parent',fig,'Style', 'text', 'Position', [w-100, h+230, 50, 20], ...
        'String', 'Width');
    txt_Y = uicontrol('Parent',fig,'Style', 'text', 'Position', [w-100, h+260, 50, 20], ...
        'String', 'Y');
    txt_X = uicontrol('Parent',fig,'Style', 'text', 'Position', [w-100, h+290, 50, 20], ...
        'String', 'X');
    txt_ftsz = uicontrol('Parent',fig,'Style', 'text', 'Position', [w-100, h+170, 50, 20], ...
        'String', 'Font size');

    label_L = uicontrol('Parent',fig,'Style', 'text', 'Position', [w+100, h+200, 50, 20], ...
        'String', num2str(L));
    label_W = uicontrol('Parent',fig,'Style', 'text', 'Position', [w+100, h+230, 50, 20], ...
        'String', num2str(W));
    label_Y = uicontrol('Parent',fig,'Style', 'text', 'Position', [w+100, h+260, 50, 20], ...
        'String', num2str(Y));
    label_X = uicontrol('Parent',fig,'Style', 'text', 'Position', [w+100, h+290, 50, 20], ...
        'String', num2str(X));
    label_ftsz = uicontrol('Parent',fig,'Style', 'text', 'Position', [w+100, h+170, 50, 20], ...
        'String', num2str(h1.FontSize));


    uicontrol('Parent',fig,'Style', 'pushbutton','String',"Start Analyzing",'Position', [w+100, h+50, 100, 30],'Callback',@SetT)
    uiwait(fig)

    function updateImage(~,~)
        L = round(get(slider_L, 'Value'));
        W = round(get(slider_W, 'Value'));
        X = round(get(slider_X, 'Value'));
        Y = round(get(slider_Y, 'Value'));
        ftsz = round(get(slider_ftsz, 'Value'));
        h1.FontSize = ftsz;
        h2.FontSize = ftsz;
        xmin = max(0,X);
        ymin = max(0,Y);
        xmax = min(X+W,size(Peak_Img,2));
        ymax = min(Y+L,size(Peak_Img,1));

        set(label_L, 'String', L);
        set(label_W, 'String', W);
        set(label_X, 'String', X);
        set(label_Y, 'String', Y);
        set(label_ftsz, 'String', ftsz);
        
        set(ax1_1, 'Xlim', [xmin,xmax]);
        set(ax1_1, 'Ylim', [ymin,ymax]);
        set(ax2_1, 'Xlim', [xmin,xmax]);
        set(ax2_1, 'Ylim', [ymin,ymax]);
        pbaspect(ax2_1,[xmax-xmin ymax-ymin 1])
        pbaspect(ax2_2,[xmax-xmin ymax-ymin 1])        
        pbaspect(ax2_3,[xmax-xmin ymax-ymin 1])        
        pbaspect(ax2_4,[xmax-xmin ymax-ymin 1])
        pbaspect(ax1_1,[xmax-xmin ymax-ymin 1])
        pbaspect(ax1_2,[xmax-xmin ymax-ymin 1])        
        pbaspect(ax1_3,[xmax-xmin ymax-ymin 1])        
        % x_vector = [xmax-Size_pixels-20 xmax-20];   % for plotting
        x_vector = [xmin+ 20 xmin+Size_pixels+20];   % for plotting
        y_vector = [ymax-20 ymax-20];   % for plotting
        L1.XData =  x_vector;
        L1.YData =  y_vector;
        L2.XData =  x_vector;
        L2.YData =  y_vector;
        T1.Position = [xmax-Size_pixels ymax-50 0];
        T2.Position = [xmax-Size_pixels ymax-50 0];
    end

    function SetT(~,~)
        delete(fig)
    end
end

function exp_dist = FitHistogram(W,param_plot,param,IC)
    %%% For Multi-exponential
    %%% Parameters are:  Amp   tau-1  sigma  ... 
    %%% Parameters are: [b(1)   b(2)  b(3)   ... 
    Costfn = @(b) sum((get_dist(b,param_plot,param) - W).^2);    
    opts = optimoptions('fmincon', 'Algorithm','interior-point','MaxFunEvals', 50000, 'MaxIter', 10000, 'StepTolerance', 1e-10,'Display', 'off');
    problem = createOptimProblem('fmincon','x0',IC,'objective',Costfn,'lb',[0,0,0,0,0,0],'ub',[10*IC(1),10,0.2,10*IC(1),10,0.2],'Aeq', [], 'beq', [],'options', opts);
    B_final = fmincon(problem);
    exp_dist = get_dist(B_final,param_plot,param);
end

function exp_dist = get_dist(B,param_plot,param)
    tau_vector = param_plot.hbins_tau;
    tau_vector(end) = [];
    exp_dist = zeros(size(tau_vector));
    for j = 1:param.order
	    amp = B(1+(j-1)*3);        
	    mu = B(2+(j-1)*3);
        sigma = B(3+(j-1)*3);
        pdf_temp = amp*normpdf(tau_vector,mu,sigma);
        exp_dist = exp_dist+pdf_temp;
   end
end

function [histw, histv] = histwv(v, w, min, max, bins)

    %Inputs: 
    %vv - values
    %ww - weights
    %minV - minimum value
    %maxV - max value
    %bins - number of bins (inclusive)
    
    %Outputs:
    %histw - weighted histogram
    %histv (optional) - histogram of values    
   
    delta = (max-min)/(bins-1);
    subs = round((v-min)/delta)+1;
    
    histw = accumarray(subs(:),w(:),[bins,1])';
    if nargout == 2
        histv = accumarray(subs(:),1,[bins,1])';
    end
end

function fig_name = get_file_name(fig_name_template,param_plot)
    fig_name_temp = strcat(fig_name_template,param_plot.figformat);
    i = 1;
    while(isfile(fig_name_temp))
        fig_name_temp = strcat(fig_name_template ,'_',num2str(i),param_plot.figformat);
        i = i+1;
    end
    fig_name = fig_name_temp;
end
