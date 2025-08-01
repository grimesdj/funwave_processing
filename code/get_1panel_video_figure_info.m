function [fig,ax0,cx0,cx1,ps,ppos0,pos] = get_1panel_video_figure_info(alims)


xm = 5;
ym = 4;
ag = 0.5;
ar = (alims(2)-alims(1))/(alims(4)-alims(3));
pw = 24;
ph = pw*ar;

ppos0 = [xm ym pw ph];
cpos0 = [xm+1*pw/8 ym+ph+ag/2 pw/5 1.25*ag]
cpos1 = [xm+3*pw/4 ym+ph+ag/2 pw/5 1.25*ag];

ps = [1.5*xm+pw, 1.3*ym+ph+3*ag];

fig = figure('units','centimeters','color','w');
pos = get(fig,'position');
pos(3:4)=ps;
set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps])

ax0 = axes('units','centimeter','position',ppos0);
cx0 = axes('units','centimeter','position',cpos0);
cx1 = axes('units','centimeter','position',cpos1);

end

