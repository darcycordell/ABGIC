function b = clean_mag_data(b)
%Data located in 04_MAG_DATA folder
%%
%[b] = load_mag_data;

%Delete duplicates
b(cellfun(@isempty,({b(:).x}))) = [];
%%
% Check for "bad" points
winlength = min([b.nt 3600]);
is = 1;tidx = 1:winlength;%length(b(is).times);
while 1
    
    plot(b(is).times(tidx),b(is).x(tidx)-mean(b(is).x(tidx),'omitnan'),'-b'); hold on
    plot(b(is).times(tidx),b(is).y(tidx)-mean(b(is).y(tidx),'omitnan'),'-r');
    ylabel('B (nT)')
    title(['Magnetic Time Series for ',upper(b(is).site)]);
    
    mmenu = menu('','Edit Points (Value)','Edit Points (Time)','Next Window','Change Window Length','Next Site','Previous Site','Delete Sites','Save','Exit');
    
    if mmenu == 1
    
    [~,y] = ginput(2);

    %if click~=3
        ax = gca;

%Option to cut by mag value (faster for removing outliers)

        val1 = num2ruler(y(1),ax.YAxis);
        val2 = num2ruler(y(2),ax.YAxis);

        id1 = find(b(is).x-mean(b(is).x,'omitnan') > min([val1, val2]) & b(is).x-mean(b(is).x,'omitnan') < max([val1, val2]));
        id2 = find(b(is).y-mean(b(is).y,'omitnan') > min([val1, val2]) & b(is).y-mean(b(is).y,'omitnan') < max([val1, val2]));

        id = unique([id1 id2]);

        b(is).x(id) = NaN;
        b(is).y(id) = NaN;
        b(is).z(id) = NaN;


        b(is).x = inpaint_nans(b(is).x,4);
        b(is).y = inpaint_nans(b(is).y,4);
        b(is).z = inpaint_nans(b(is).z,4);
        close(gcf);

    elseif mmenu == 2

        [x,~] = ginput(2);


        ax = gca;
        %Option to cut by time
        t1 = num2ruler(x(1),ax.XAxis);
        t2 = num2ruler(x(2),ax.XAxis);
        [~ ,id1] = min(abs(datenum(t1)-datenum(b(is).times)));
        [~ ,id2] = min(abs(datenum(t2)-datenum(b(is).times)));
% 
%         b(is).x(id1:id2) = linspace(b(is).x(id1-1),b(is).x(id2+1),length(id1:id2));
%         b(is).y(id1:id2) = linspace(b(is).y(id1-1),b(is).y(id2+1),length(id1:id2));
%         b(is).z(id1:id2) = linspace(b(is).z(id1-1),b(is).z(id2+1),length(id1:id2));

        b(is).x(id1:id2) = NaN;
        b(is).y(id1:id2) = NaN;
        b(is).z(id1:id2) = NaN;

        b(is).x = inpaint_nans(b(is).x,4);
        b(is).y = inpaint_nans(b(is).y,4);
        b(is).z = inpaint_nans(b(is).z,4);
        close(gcf);


    elseif mmenu == 3
        tidx = tidx+winlength;
        
        if tidx(end)>b(is).nt
            tidx = b(is).nt-winlength:b(is).nt;
        end
        
        if length(tidx)>b(is).nt
            tidx = 1:b(is).nt;
        end
        
        close(gcf);
        
    elseif mmenu == 4
        
        prompt = {'New Window Length'};
        titles  = '';
        def = {'3600'};
        tmp = inputdlg(prompt,titles,1,def); 
        winlength = str2double(tmp{1});
        
        tidx = tidx(1):tidx(1)+winlength;
        
        if tidx(end)>b(is).nt
           tidx = tidx(1):b(is).nt;
        end
        
        close(gcf);
        
    elseif mmenu == 5 %next site
        is = is+1;
        
        tidx = 1:winlength;
        
        if is>length(b)
            is = length(b);
        end
        
        if tidx(end)>b(is).nt
           tidx = tidx(1):b(is).nt;
        end
        
        close(gcf);
        
    elseif mmenu == 6 %previous site

        is = is-1;
        
        tidx = 1:winlength;
        
        if is<1
            is = 1;
        end
        
        if tidx(end)>b(is).nt
           tidx = tidx(1):b(is).nt;
        end
        
        close(gcf);

    elseif mmenu == 7 %delete site

        table({b.site});
    
        [sel,val] = listdlg('PromptString',['Select stations to delete (',num2str(length(b)),' total) '],'ListString',{b.site},'Name','Delete Stations','ListSize',[300 300]);
   
        if val

            b(sel) = [];

        end
        close(gcf)

        
    elseif mmenu == 8 %save
        t0 = b(1).times(1);
        t1 = datetime;
        savestr = [num2str(year(t0)) num2str(month(t0),'%02.f') num2str(day(t0),'%02.f') '_edit_v_',...
            num2str(year(t1)) num2str(month(t1),'%02.f') num2str(day(t1),'%02.f') '.mat'];
        uisave({'b'},savestr);


    else
        break
    end


end


    

