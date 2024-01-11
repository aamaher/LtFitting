count = 0;
hbins_tau = 0.1:0.1:15;
E = hbins_tau;

%% Getting N_Obj, comment out later
% N_Obj = 1;
% for idx_Day = 1:N_days    
%     Mask_temp = Mask{idx_Day,idx_Measurement};
%     N_Obj_temp = max(Label_Mask_full(:));
%     if (N_Obj_temp > N_Obj)
%         N_Obj = N_Obj_temp;
%     end
% end

%%  Plot for each label separately

for idx_Obj = 1:N_Obj


    FullMask_temp = Label_Mask_full;
    FullMask_temp ((Label_Mask_full ~= idx_Obj) & (Label_Mask_full ~= 0)) = 1;
    FullMask_temp (Label_Mask_full == idx_Obj) = 2;
    if (sum(FullMask_temp == 2) == 0)
          disp(strcat("Not enough points for label: ",num2str(idx_Obj)));
          continue;
    end
    meanvalue = zeros(N_days,1);
    idx_mask = 1;
    Points_mask = 0;
    while (~Points_mask)
        Mask_Temp = Mask{idx_mask};
        Mask_Temp(Mask_Temp ~= idx_Obj) = 0;
        Mask_Temp(Mask_Temp == idx_Obj) = 1;
        Points_mask = sum(Mask_Temp(:));
        if (~Points_mask)
            if (idx_mask <= N_days)
                idx_mask = idx_mask + 1;
            else
                disp(strcat("Weird error in label: ",num2str(idx_Obj)));
                continue;
            end
        end
    end
    stats = regionprops(Mask_Temp, 'Centroid');
    center = stats.Centroid;        %%% Center of the Mask for the current label

    tau_Obj = tau_map_all{idx_mask};
    tau_Obj(~Mask_Temp) = NaN;
    gcf = figure(idx_Obj*(idx_Measurement-1));
    gcf.Position = [50 50 1000 600];
    ax1 = axes(gcf,'Position',[0.05 0.22 0.6 0.6]);
    ax2 = copyobj(ax1,gcf);
    imagesc(ax1,Img_Big_Full); %% Plot the image
    imagesc(ax2,FullMask_temp, 'AlphaData', 0.5*FullMask_temp);

    colormap(ax1,'gray')
    colormap(ax2,gray(3))
    set(ax1,'XColor', 'none','YColor','none')
    set(ax2,'XColor', 'none','YColor','none')

    ax2.UserData = linkprop([ax1,ax2],{'Position','InnerPosition','DataAspectRatio','xtick','ytick','ydir','xdir','xlim','ylim'});
    ax2.Visible = 'off';
    title(['Object Number:' num2str(idx_Obj) ', at day:' num2str(1)],'Fontsize',20);
    xmin = find(sum(Mask_Temp,1), 1, 'first');
    ymin = find(sum(Mask_Temp,2), 1, 'first');
    xmax = find(sum(Mask_Temp,1), 1, 'last');
    ymax = find(sum(Mask_Temp,2), 1, 'last');
    dx = xmax - xmin;
    dy = ymax - ymin;
    % rectangle('Position', [center(1),center(2), 120, 120],'EdgeColor', 'r', 'LineWidth', 2);
    
    %%% Zoomed in figure
    ax4 = axes('Position',[0.67 0.2 0.28 0.28]);hold on;
    ax5 = copyobj(ax4,gcf);
    imagesc(ax4,Peak_Img_all{idx_mask}); %% Plot the image
    imagesc(ax5,tau_Obj, 'AlphaData', Mask_Temp);
    xlim([max(center(1)-dx,0),min(center(1)+dx,size(Peak_Img_all{1},2))])
    ylim([max(center(2)-dy,0),min(center(2)+dy,size(Peak_Img_all{1},1))])
    colormap(ax4,'gray')
    colormap(ax5,jet(100))


    [V,E] = histcounts(tau_Obj,param_plot.hbins_tau);
    W = V/sum(V);
    dE = diff(E(1:2));
    tau_mean = sum((E(1:end-1)+dE/2).*W);
    tau_std = sqrt(sum(W.*(E(1:end-1) + dE/2 - tau_mean).^2));


    clim(ax5,[max(0,tau_mean-3*tau_std),tau_mean+3*tau_std]);

    c = colorbar(ax5);
    c.Label.String = 'Lifetime (ns)';
    c.Label.FontSize = 14;

    set(ax4,'XColor', 'none','YColor','none')
    set(ax5,'XColor', 'none','YColor','none')
    ax5.UserData = linkprop([ax4,ax5],{'Position','InnerPosition','DataAspectRatio','xtick','ytick','ydir','xdir','xlim','ylim'});
    ax5.Visible = 'off';
    title(strcat('Day ',num2str(idx_mask)))
    %%% Histogram
    ax3 = axes('Position',[0.72 0.6 0.2 0.28]);hold on;
    leg_str  = cell(N_days,1);
    for idx_day = 1:N_days

        Mask_Temp = Mask{idx_day};
        Mask_Temp(Mask_Temp ~= idx_Obj) = 0;
        Mask_Temp(Mask_Temp == idx_Obj) = 1;

        tau_Obj = tau_map_all{idx_day};
        tau_Obj(~Mask_Temp) = NaN;

        [V,E] = histcounts(tau_Obj,hbins_tau);
        W = V/sum(V);
        meanvalue(idx_day) = sum(E(1:end-1).*W);

        plot(ax3,E(1:end-1),W,'LineWidth',2);
        set(ax3,'XLim',[0,5])
        hx = xlabel('Lifetime (ns)');
        hy = ylabel('Distribution');
        % xline(meanvalue,'LineWidth',1.5)
        leg_str{idx_day} = strcat('Day ',num2str(idx_day));
        if (idx_day == 1)
            title_str = strcat("Mean: ",num2str(meanvalue(idx_day)));
        else
            title_str = strcat(title_str,"->",num2str(meanvalue(idx_day)));
        end
        
    end
    leg = legend (string(leg_str),'location','northwest');
    title(title_str)
end

