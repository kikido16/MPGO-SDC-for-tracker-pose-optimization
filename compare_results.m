init_trans=[1.80807630477183,5.41040945208734,8.99389831015270];

cov_trans=[1.58964433782277,4.73906204230771,7.96973938767133];
rf_trans=[1.57282289712453,4.72342341798303,7.89663217520124];

% cov_trans=[1.19833547469887,3.59862511336060,6.20153503991649];
% rf_trans=[1.15235390762402,3.76289006938175,6.52032444985688];


ref_trans=[1.70633198377182	4.99609838780641	8.30349639321889];

init_rot=[0.0194320997700625,0.0582945222226215,0.0971516097892983];

cov_rot=[0.0125389052888270,0.0370398743463621,0.0647918546328938];

rf_rot=[0.0121622094106117,0.0365298115730793,0.0615214959964078];
% cov_rot=[0.0119458918417615,0.0358054789240668,0.0596988660574912];
% 
% rf_rot=[0.0119774067176170,0.0360997734336477,0.0627918138127060];

ref_rot=[0.0256864741311504	0.0764088236872382	0.127207701526040];

init_trans2=[4.75836917897213,14.3142048931052,23.9170332393284];

cov_trans2=[2.52934768149586,7.46818955492227,12.6030117100230];

rf_trans2=[2.59432131247662,7.88864562836939,13.0145341812938];
% 
% cov_trans2=[2.11436203800630,7.13721778490923,13.3871386813439];
% 
% rf_trans2=[3.02835298907545,7.84469735776413,16.6898106175909];

ref_trans2=[3.74061033587412	11.0675342607997	18.4365590034041];

init_rot2=[0.0309077287818817,0.0928401819366730,0.154911381732232];

cov_rot2=[0.0175444134516178,0.0520109281600747,0.0864524895110184];

rf_rot2=[0.0193295222063513,0.0584497992355087,0.0978915479085247];
% cov_rot2=[0.0174613345794775,0.0515325424181139,0.0888525062164203];
% 
% rf_rot2=[0.0235943222340381,0.0637719674400301,0.124021006189973];


ref_rot2=[0.0316802297323655	0.0942872894889159	0.157086348429735];



%%绘制3D柱状图
trans_result_data=[init_trans',ref_trans',cov_trans',rf_trans'];

figure(1)
b1=bar3(trans_result_data,0.7);
hXLabel = xlabel('Noise level groups','FontWeight','bold','FontSize',15,'Rotation',8,'Position',[2 8 3]);
hYLabel = ylabel('Method','FontWeight','bold','FontSize',15,'Rotation',-30,'Position',[-0.15 3.5 -3.5]);
hXLabel = xlabel('Method','FontWeight','bold','FontSize',15,'Rotation',8,'Position',[2.65 -0.08 -1.1]);
hYLabel = ylabel('Noise level groups','FontWeight','bold','FontSize',15,'Rotation',-15,'Position',[4.1 1.65 -1.6]);
hZLabel = zlabel('Average Et (mm)','FontWeight','bold','FontSize',15);
for k = 1:length(b1)
    zdata = b1(k).ZData;
    b1(k).CData = zdata;
    b1(k).FaceColor = 'interp';
end
% colormap(slanCM('rainbow-iso'));
colormap(linspecer);
set(gca, 'Box', 'on', ...                                                         % 边框
    'LineWidth', 1, 'GridLineStyle', '-',...                                   % 坐标轴线宽
    'XGrid', 'on', 'YGrid', 'on','ZGrid', 'on', ...                          % 网格
    'TickDir', 'out', 'TickLength', [.015 .015], ...                           % 刻度
    'XMinorTick', 'off', 'YMinorTick', 'off',  'ZMinorTick', 'off',...         % 小刻度
    'XColor', [.1 .1 .1],  'YColor',[.1 .1 .1], 'ZColor', [.1 .1 .1],...       % 坐标轴颜色
    'Yticklabel',{'Low','Medium','High'},...                                                    % X坐标轴刻度标签
    'Xticklabel',{'BF','SB','Cov','RF'},'FontSize',20);
view(-140,10);

rot_result_data=[init_rot',ref_rot',cov_rot',rf_rot'];

figure(2)
b2=bar3(rot_result_data,0.7);
% hXLabel = xlabel('Method','FontWeight','bold','FontSize',15,'Rotation',8,'Position',[2.6 0 -0.02]);
% hYLabel = ylabel('Noise level groups','FontWeight','bold','FontSize',15,'Rotation',-15,'Position',[4.1 1.65 -0.025]);
% hZLabel = zlabel('Average Er (\circ)','FontWeight','bold','FontSize',15);
for k = 1:length(b2)
    zdata = b2(k).ZData;
    b2(k).CData = zdata;
    b2(k).FaceColor = 'interp';
end
% colormap(slanCM('rainbow-iso'));
colormap(linspecer);
set(gca, 'Box', 'on', ...                                                         % 边框
    'LineWidth', 1, 'GridLineStyle', '-',...                                   % 坐标轴线宽
    'XGrid', 'on', 'YGrid', 'on','ZGrid', 'on', ...                          % 网格
    'TickDir', 'out', 'TickLength', [.015 .015], ...                           % 刻度
    'XMinorTick', 'off', 'YMinorTick', 'off',  'ZMinorTick', 'off',...         % 小刻度
    'XColor', [.1 .1 .1],  'YColor',[.1 .1 .1], 'ZColor', [.1 .1 .1],...       % 坐标轴颜色
    'Yticklabel',{'Low','Medium','High'},...                                                    % X坐标轴刻度标签
    'Xticklabel',{'BF','SB','Cov','RF'},'FontSize',20);
view(-140,10);

trans2_result_data=[init_trans2',ref_trans2',cov_trans2',rf_trans2'];

figure(3)
b3=bar3(trans2_result_data,0.7);
% hXLabel = xlabel('Method','FontWeight','bold','FontSize',15,'Rotation',8,'Position',[2.65 -0.08 -2.75]);
% hYLabel = ylabel('Noise level groups','FontWeight','bold','FontSize',15,'Rotation',-15,'Position',[4.1 1.65 -4]);
% hZLabel = zlabel('Average Et (mm)','FontWeight','bold','FontSize',15);
for k = 1:length(b3)
    zdata = b3(k).ZData;
    b3(k).CData = zdata;
    b3(k).FaceColor = 'interp';
end
% colormap(slanCM('rainbow-iso'));
colormap(linspecer);
set(gca, 'Box', 'on', ...                                                         % 边框
    'LineWidth', 1, 'GridLineStyle', '-',...                                   % 坐标轴线宽
    'XGrid', 'on', 'YGrid', 'on','ZGrid', 'on', ...                          % 网格
    'TickDir', 'out', 'TickLength', [.015 .015], ...                           % 刻度
    'XMinorTick', 'off', 'YMinorTick', 'off',  'ZMinorTick', 'off',...         % 小刻度
    'XColor', [.1 .1 .1],  'YColor',[.1 .1 .1], 'ZColor', [.1 .1 .1],...       % 坐标轴颜色
    'Yticklabel',{'Low','Medium','High'},...                                                    % X坐标轴刻度标签
    'Xticklabel',{'BF','SB','Cov','RF'},'FontSize',20);
view(-140,10);

rot2_result_data=[init_rot2',ref_rot2',cov_rot2',rf_rot2'];

figure(4)
b4=bar3(rot2_result_data,0.7);
% hXLabel = xlabel('Method','FontWeight','bold','FontSize',15,'Rotation',8,'Position',[2.6 -0.006 -0.024]);
% hYLabel = ylabel('Noise level groups','FontWeight','bold','FontSize',15,'Rotation',-15,'Position',[4.1 1.65 -0.032]);
% hZLabel = zlabel('Average Er (\circ)','FontWeight','bold','FontSize',15);
for k = 1:length(b4)
    zdata = b4(k).ZData;
    b4(k).CData = zdata;
    b4(k).FaceColor = 'interp';
end
% colormap(slanCM('rainbow-iso'));
colormap(linspecer);
set(gca, 'Box', 'on', ...                                                         % 边框
    'LineWidth', 1, 'GridLineStyle', '-',...                                   % 坐标轴线宽
    'XGrid', 'on', 'YGrid', 'on','ZGrid', 'on', ...                          % 网格
    'TickDir', 'out', 'TickLength', [.015 .015], ...                           % 刻度
    'XMinorTick', 'off', 'YMinorTick', 'off',  'ZMinorTick', 'off',...         % 小刻度
    'XColor', [.1 .1 .1],  'YColor',[.1 .1 .1], 'ZColor', [.1 .1 .1],...       % 坐标轴颜色
    'Yticklabel',{'Low','Medium','High'},...                                                    % X坐标轴刻度标签
    'Xticklabel',{'BF','SB','Cov','RF'},'FontSize',20);
view(-140,10);
