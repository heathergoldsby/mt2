%% Initialization
addpath /Users/dk/research/prj/distributed-control/doc/gecco2010/math -begin
cd /Users/dk/research/prj/distributed-control/var

%%
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"001-structure/n_.*gen-Genchamp-AvgFit"', 'N';... % 15
    '"001-structure/ini_.*gen-Genchamp-AvgFit"', 'H1';... % 1
    '"001-structure/mhn_.*gen-Genchamp-AvgFit"', 'M1';... % 11
    },...
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('Initial');
quick_export('../doc/gecco2010/figures/001-initial.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%%
% Ok, this box plot is showing us the overall fitness ratios of all the
% different structures that we tested.  As a result of the Kruskal-Wallis
% test, which shows us the multiple-comparison of all of these approaches,
% we selected Neat, 3D-Hyperneat, and multiagent hyperneat to further test.
%
% The approaches we tested were not significantly different from NEAT (or
% each other) at the p<0.01 level.
clear d1 m;
d1 = quick_boxplot('data(:,end,2)',...
    {...
    '"001-structure/ini_.*gen-Genchamp-AvgFit"', 'H1';... % 1
    '"001-structure/ini-h_.*gen-Genchamp-AvgFit"', 'H1+';...
    '"001-structure/geo1_.*gen-Genchamp-AvgFit"', 'H2';...
    '"001-structure/geo1-h_.*gen-Genchamp-AvgFit"', 'H2+';...
    '"001-structure/geo2_.*gen-Genchamp-AvgFit"', 'H3';... % 5
    '"001-structure/geo2-h_.*gen-Genchamp-AvgFit"', 'H3+';...
    '"001-structure/spread_.*gen-Genchamp-AvgFit"', 'H4';...
    '"001-structure/3d_.*gen-Genchamp-AvgFit"', 'D1';...
    '"001-structure/3d-h_.*gen-Genchamp-AvgFit"', 'D1+';...
    '"001-structure/3d2_.*gen-Genchamp-AvgFit"', 'D2';...
    '"001-structure/3d2-h_.*gen-Genchamp-AvgFit"', 'D2+';... % 10
    '"001-structure/3d3_.*gen-Genchamp-AvgFit"', 'D3';...
    '"001-structure/3d3-h_.*gen-Genchamp-AvgFit"', 'D3+';...
    '"001-structure/3d4_.*gen-Genchamp-AvgFit"', 'D4';...
    '"001-structure/mhn_.*gen-Genchamp-AvgFit"', 'M1';...
    '"001-structure/mhn-h_.*gen-Genchamp-AvgFit"', 'M1+';... % 15
    '"001-structure/mhn3_.*gen-Genchamp-AvgFit"', 'M2';...
    '"001-structure/mhn3-h_.*gen-Genchamp-AvgFit"', 'M2+';...
    '"001-structure/rmhn_.*gen-Genchamp-AvgFit"', 'R1';...
    '"001-structure/rmhn-h_.*gen-Genchamp-AvgFit"', 'R1+';...
    '"001-structure/rmhn2_.*gen-Genchamp-AvgFit"', 'R2';... %20
    '"001-structure/rmhn2-h_.*gen-Genchamp-AvgFit"', 'R2+';...
    '"001-structure/n_.*gen-Genchamp-AvgFit"', 'N';...
    '"002-factors/sense6-n_.*gen-Genchamp-AvgFit"', 'NC';... % neat w/ cell sensors
    '"001-structure/nearest_.*gen-Genchamp-AvgFit"', 'NE';...
    '"001-structure/nohidden_.*gen-Genchamp-AvgFit"', 'NH';... %25
    '"001-structure/twoout_.*gen-Genchamp-AvgFit"', 'T';...
    '"001-structure/bias-d2_.*gen-Genchamp-AvgFit"', 'B1';...
    '"001-structure/bias-m1_.*gen-Genchamp-AvgFit"', 'B2';...
    '"001-structure/faq-h1_.*gen-Genchamp-AvgFit"', 'F1';... % 30
    '"001-structure/faq-m1_.*gen-Genchamp-AvgFit"', 'F2';...
    },...
    @(i,x) (x./1600));

xlabel('Approach')
ylabel('Fraction of max. fitness');
title('Varied neuroevolutionary approaches');
quick_export('../doc/gecco2010/figures/001-structure-box.eps', 12, 3.75);

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
m(1:30,4) = d1{4}.data(:,end,2);
m(1:30,5) = d1{5}.data(:,end,2);
m(1:30,6) = d1{6}.data(:,end,2);
m(1:30,7) = d1{7}.data(:,end,2);
m(1:30,8) = d1{8}.data(:,end,2);
m(1:30,9) = d1{9}.data(:,end,2);
m(1:30,10) = d1{10}.data(:,end,2);
m(1:30,11) = d1{11}.data(:,end,2);
m(1:30,12) = d1{12}.data(:,end,2);
m(1:30,13) = d1{13}.data(:,end,2);
m(1:30,14) = d1{14}.data(:,end,2);
m(1:30,15) = d1{15}.data(:,end,2);
m(1:30,16) = d1{16}.data(:,end,2);
m(1:30,17) = d1{17}.data(:,end,2);
m(1:30,18) = d1{18}.data(:,end,2);
m(1:30,19) = d1{19}.data(:,end,2);
m(1:30,20) = d1{20}.data(:,end,2);
m(1:30,21) = d1{21}.data(:,end,2);
m(1:30,22) = d1{22}.data(:,end,2);
m(1:30,23) = d1{23}.data(:,end,2);
m(1:30,24) = d1{24}.data(:,end,2);
m(1:30,25) = d1{25}.data(:,end,2);
m(1:30,26) = d1{26}.data(:,end,2);
m(1:30,27) = d1{27}.data(:,end,2);
m(1:30,28) = d1{28}.data(:,end,2);
m(1:30,29) = d1{29}.data(:,end,2);
m(1:30,30) = d1{30}.data(:,end,2);
m(1:30,31) = d1{31}.data(:,end,2);

[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);


%%
% Here we see that NEAT is statistically better than either of the two
% HyperNEAT approaches, while there is no statistical difference betweeen
% HyperNEAT approaches (p<0.01).
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"002-factors/n-id_.*gen-Genchamp-AvgFit"', 'N';...
    '"002-factors/3d-id_.*gen-Genchamp-AvgFit"', 'D1';...
    '"002-factors/mhn-id_.*gen-Genchamp-AvgFit"', 'M1'},...
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('+ID inputs');
quick_export('../doc/gecco2010/figures/002-id.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%%
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"002-factors/sense5-hn_.*gen-Genchamp-AvgFit"', '5-HN';...
    '"002-factors/sense6-hn_.*gen-Genchamp-AvgFit"', '6-HN';...
    '"002-factors/sense5-n_.*gen-Genchamp-AvgFit"', '5-N';...
    '"002-factors/sense6-n_.*gen-Genchamp-AvgFit"', '6-N';...
    '"001-structure/n_.*gen-Genchamp-AvgFit"', 'N';...
    },...
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('Sensor configurations');
%quick_export('../doc/gecco2010/figures/002-sense.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
m(1:30,4) = d1{2}.data(:,end,2);
m(1:30,5) = d1{3}.data(:,end,2);

[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%%
% As above, NEAT is statistically significantly better, and the two
% HyperNEAT approaches are not statistically different from each other at
% the p<0.01 level.
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"002-factors/n-rxd_.*gen-Genchamp-AvgFit"', 'N';...
    '"002-factors/3d-rxd_.*gen-Genchamp-AvgFit"', 'D1';...
    '"002-factors/mhn-rxd_.*gen-Genchamp-AvgFit"', 'M1'},...
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('+Received data');
quick_export('../doc/gecco2010/figures/002-rxd.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:29,3) = d1{3}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%%
% Same same.  NEAT is statistically significantly better, and the two
% HyperNEAT approaches are statistically different from each other at the
% p<0.01 level.
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"002-factors/n-rxnull_.*gen-Genchamp-AvgFit"', 'N';...
    '"002-factors/3d-rxnull_.*gen-Genchamp-AvgFit"', 'D1';...
    '"002-factors/mhn-rxnull_.*gen-Genchamp-AvgFit"', 'M1'},...
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('-Rx sensors');
quick_export('../doc/gecco2010/figures/002-rxnull.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%%
% NEAT is statistically significantly better without location sensors, but
% here the 3D treatment is also significantly better than the multi-agent
% treatment (at the p<0.01 level).  Interesting.
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"002-factors/n-locnull_.*gen-Genchamp-AvgFit"', 'N';...
    '"002-factors/3d-locnull_.*gen-Genchamp-AvgFit"', 'D1';...
    '"002-factors/mhn-locnull_.*gen-Genchamp-AvgFit"', 'M1'},...
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('-Loc sensors');
quick_export('../doc/gecco2010/figures/002-locnull.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%%
% This is showing us that yes, the sheer number of cells that are available
% to be occupied influences fitness.  The difference between 6 & 8 is not
% significant, however the difference between 6 & 10 is (at the p<0.01
% level).
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"001-structure/3d2_.*gen-Genchamp-AvgFit"', '6';...
    '"003-gridsize/3d-8_.*gen-Genchamp-AvgFit"', '8';...
    '"003-gridsize/3d-10_.*gen-Genchamp-AvgFit"', '10'},...
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('Varying grid size, D1');
quick_export('../doc/gecco2010/figures/003-grid-3d.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%% hjg - 003- neat - vary grid size
% Unfortunately, all results here are significantly different from each
% other.  This means that as we increase the size of the grid, networks are
% more likely to occupy different cells via chance alone.
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"001-structure/n_.*gen-Genchamp-AvgFit"', '6';...
    '"003-gridsize/n-8_.*gen-Genchamp-AvgFit"', '8';...
    '"003-gridsize/n-10_.*gen-Genchamp-AvgFit"', '10'},...
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('Varying grid size, NEAT');
quick_export('../doc/gecco2010/figures/003-grid-n.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%% hjg - 003 - does radio range affect performance?
% In short, yes.  Radios of range 6 & 8 are not significantly different
% from each other, however {6,8} are significantly different from 10, with
% longer-range radios performing better; this makes sense, as it's easier
% to maintain a connected network.
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"001-structure/3d2_.*gen-Genchamp-AvgFit"', '10';...
    '"003-gridsize/3d-rad-8_.*gen-Genchamp-AvgFit"', '8';...
    '"003-gridsize/3d-rad-6_.*gen-Genchamp-AvgFit"', '6'},...
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('Varying radio range, D1');
quick_export('../doc/gecco2010/figures/003-rad-3d.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%% hjg - 003 - does radio range affect performance?
% Yes, but not as strongly.  6 and 10 are significantly different from each
% other, but neither 6 & 8 nor 8 & 10 are.
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"001-structure/n_.*gen-Genchamp-AvgFit"', '10';...
    '"003-gridsize/n-rad-8_.*gen-Genchamp-AvgFit"', '8';...
    '"003-gridsize/n-rad-6_.*gen-Genchamp-AvgFit"', '6'},...
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('Varying radio range, NEAT');
quick_export('../doc/gecco2010/figures/003-rad-n.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%% hjg - 003 - neat contingent on max velo
% Significant differences between 8 and 4, with higher velocities being
% better.
%
% In general, when the 003 series of experiments have shown us is that both
% hyperneat and neat take advantage of similar features of the
% environment... perhaps we can look at the ratios of fitness improvements
% for each of these systems?
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"001-structure/n_.*gen-Genchamp-AvgFit"', '8';...
    '"003-gridsize/n-low-6_.*gen-Genchamp-AvgFit"', '6';...
    '"003-gridsize/n-low-4_.*gen-Genchamp-AvgFit"', '4'},...
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('Varying maximum velocity, NEAT');
quick_export('../doc/gecco2010/figures/003-velocity-n.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%%
% As expected, what we're seeing here is that increasing the size of the
% network increases the difficulty of the problem.
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"004-netsize/n2_.*gen-Genchamp-AvgFit"', '2';...
    '"004-netsize/n4_.*gen-Genchamp-AvgFit"', '4';...
    '"004-netsize/n8_.*gen-Genchamp-AvgFit"', '8';...    
    '"001-structure/n_.*gen-Genchamp-AvgFit"', '16';...
    '"004-netsize/n32_.*gen-Genchamp-AvgFit"', '32'},...
    @(i,x) (x./(2^i*100)));
axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('Varying network size, NEAT');
quick_export('../doc/gecco2010/figures/004-size-n.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
m(1:30,4) = d1{4}.data(:,end,2);
m(1:30,5) = d1{5}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%%
x = [mean(d1{1}.data(:,end,2))/200;...
    mean(d1{2}.data(:,end,2))/400;...
    mean(d1{3}.data(:,end,2))/800;...
    mean(d1{4}.data(:,end,2))/1600;...
    mean(d1{5}.data(:,end,2))/3200];

y = zeros(size(x,1)-1,1);
for i=1:size(x,1)-1
    lx = mean(x(i)-x(i+1));
    y(i) = lx;
end
display(mean(y));

%%
% As expected, what we're seeing here is that increasing the size of the
% network increases the difficulty of the problem.
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"004-netsize/3d2_.*gen-Genchamp-AvgFit"', '2';...
    '"004-netsize/3d4_.*gen-Genchamp-AvgFit"', '4';...
    '"004-netsize/3d8_.*gen-Genchamp-AvgFit"', '8';...    
    '"001-structure/3d2_.*gen-Genchamp-AvgFit"', '16';...
    '"004-netsize/3d32_.*gen-Genchamp-AvgFit"', '32'},...
    @(i,x) (x./(2^i*100)));
axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('Varying network size, 3DHN');
quick_export('../doc/gecco2010/figures/004-size-3d.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
m(1:30,4) = d1{4}.data(:,end,2);
m(1:29,5) = d1{5}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%%
% Again, increasing the size of the network increases the difficulty of the
% problem.  NEAT out-performed multiagent hyperneat here, too.
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"004-netsize/mhn2_.*gen-Genchamp-AvgFit"', '2';...
    '"004-netsize/mhn4_.*gen-Genchamp-AvgFit"', '4';...
    '"004-netsize/mhn8_.*gen-Genchamp-AvgFit"', '8';...    
    '"001-structure/mhn_.*gen-Genchamp-AvgFit"', '16';...
    '"004-netsize/mhn32_.*gen-Genchamp-AvgFit"', '32'},...
    @(i,x) (x./(2^i*100)));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('Varying network size, Multi-agent HyperNEAT');
quick_export('../doc/gecco2010/figures/004-size-mhn.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
m(1:30,4) = d1{4}.data(:,end,2);
m(1:30,5) = d1{5}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%%
clear d1 m;
d1 = quick_multiload(...
    {...
    '"001-structure/n_.*gen-Genchamp-AvgFit"', '16';...
    '"004-netsize/n32_.*gen-Genchamp-AvgFit"', '32';...
    '"001-structure/3d2_.*gen-Genchamp-AvgFit"', '16';...
    '"004-netsize/3d32_.*gen-Genchamp-AvgFit"', '32';...
    '"001-structure/mhn_.*gen-Genchamp-AvgFit"', '16';...
    '"004-netsize/mhn32_.*gen-Genchamp-AvgFit"', '32';...
    });
m(1:29,1) = d1{1}.data(1:29,end,2) - d1{2}.data(1:29,end,2);
m(1:29,2) = d1{3}.data(1:29,end,2) - d1{4}.data(1:29,end,2);
m(1:29,3) = d1{5}.data(1:29,end,2) - d1{6}.data(1:29,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);



%%
% Here we're increasing the number of radio sensors, and we find no
% significant differences between corresponding treatments (Neat is not
% sig. diff. from Neat+2, etc).
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"001-structure/n_.*gen-Genchamp-AvgFit"', 'N';...
    '"001-structure/3d2_.*gen-Genchamp-AvgFit"', 'D1';...
    '"001-structure/mhn_.*gen-Genchamp-AvgFit"', 'M1';...
    '"005-hifi/n6_.*gen-Genchamp-AvgFit"', 'N+2';...
    '"005-hifi/3d6_.*gen-Genchamp-AvgFit"', 'D1+2';...
    '"005-hifi/mhn6_.*gen-Genchamp-AvgFit"', 'M1+2'},...    
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('+2 sensors');
quick_export('../doc/gecco2010/figures/005-hifi.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
m(1:30,4) = d1{4}.data(:,end,2);
m(1:30,5) = d1{5}.data(:,end,2);
m(1:30,6) = d1{6}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%%
% As above, here we're increasing the number of radio sensors, but also
% changing the meaning of the sensors to be the number of nodes in that
% 'pie slice' region.  In most cases, this is approach is significantly
% worse than using signal strength.
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"001-structure/n_.*gen-Genchamp-AvgFit"', 'N';...
    '"001-structure/3d2_.*gen-Genchamp-AvgFit"', 'D1';...
    '"001-structure/mhn_.*gen-Genchamp-AvgFit"', 'M1';...
    '"005-hifi/n6b_.*gen-Genchamp-AvgFit"', 'N+2';...
    '"005-hifi/3d6b_.*gen-Genchamp-AvgFit"', 'D1+2';...
    '"005-hifi/mhn6b_.*gen-Genchamp-AvgFit"', 'M1+2'},...    
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('+2 sensors, using number of nodes instead of signal strength');
quick_export('../doc/gecco2010/figures/005-hifi-b.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
m(1:30,4) = d1{4}.data(:,end,2);
m(1:30,5) = d1{5}.data(:,end,2);
m(1:30,6) = d1{6}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%% 
% Strangely enough, these plots show that the original normalization was
% good.  At this point, I have no idea what NEAT/HyperNEAT are doing...
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"001-structure/n_.*gen-Genchamp-AvgFit"', 'N';...
    '"006-renorm/n_.*gen-Genchamp-AvgFit"', 'N-renorm'},...    
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('NEAT, renormalization');
%quick_export('../doc/gecco2010/figures/005-hifi-b.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"001-structure/3d2_.*gen-Genchamp-AvgFit"', '3D';...
    '"006-renorm/3d2_.*gen-Genchamp-AvgFit"', '3D-renorm'},...    
    @(i,x) (x./1600));

axis([0 200 0 1.05]);
xlabel('Generations')
ylabel('Fraction of max. fitness');
title('3D, renormalization');
%quick_export('../doc/gecco2010/figures/005-hifi-b.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%%
% Ok, weird.  This is on the circular coverage problem, and it's showing us
% that the maximum fitness of Multi-agent HyperNEAT is greater than NEAT.
% However, these fitnesses are NOT significantly different at the p<0.01
% level.  Strategies are vastly different.
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"007-fit/circ-n_.*gen-Genchamp-AvgFit"', 'N';...
    '"007-fit/circ-hn3_.*gen-Genchamp-AvgFit"', 'HN3';...    
    '"007-fit/circ-mhn_.*gen-Genchamp-AvgFit"', 'M1'});

axis([0 200 0 1.05]);
axis 'auto y';
xlabel('Generations')
ylabel('Fitness');
title('Coverage of circular region');
quick_export('../doc/gecco2010/figures/007-circ.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%%
% This is broken??
clear d1 m;
d1 = quick_plot('data(1,:,1)', 'data(:,:,2)',...
    {...
    '"007-fit/cc-n_.*gen-Genchamp-AvgFit"', 'N';...
    '"007-fit/cc-hn3_.*gen-Genchamp-AvgFit"', 'HN3';...    
    '"007-fit/cc-mhn_.*gen-Genchamp-AvgFit"', 'M1'});

axis([0 200 0 1.05]);
axis 'auto y';
xlabel('Generations')
ylabel('Fitness');
title('Coverage of circular region');
%quick_export('../doc/gecco2010/figures/007-cc.eps');

m(1:30,1) = d1{1}.data(:,end,2);
m(1:30,2) = d1{2}.data(:,end,2);
m(1:30,3) = d1{3}.data(:,end,2);
[p, table, stats] = kruskalwallis(m);
multcompare(stats,'alpha', 0.01);

%%
cd /Users/dk/research/prj/distributed-control/var
d1 = quick_multiload({...
    '"001-structure/ini_.*gen-Genchamp-AvgFit"', 'HN1';... % 1
    '"001-structure/ini-h_.*gen-Genchamp-AvgFit"', 'HN1+H';...
    '"001-structure/geo1_.*gen-Genchamp-AvgFit"', 'HN2';...
    '"001-structure/geo1-h_.*gen-Genchamp-AvgFit"', 'HN2+H';...
    '"001-structure/geo2_.*gen-Genchamp-AvgFit"', 'HN3';... % 5
    '"001-structure/geo2-h_.*gen-Genchamp-AvgFit"', 'HN3';...
    '"001-structure/3d_.*gen-Genchamp-AvgFit"', '3D1';...
    '"001-structure/3d-h_.*gen-Genchamp-AvgFit"', '3D1+H';...
    '"001-structure/3d2_.*gen-Genchamp-AvgFit"', '3D2';...
    '"001-structure/3d2-h_.*gen-Genchamp-AvgFit"', '3D2+H';... % 10
    '"001-structure/3d3_.*gen-Genchamp-AvgFit"', '3D3';...
    '"001-structure/3d3-h_.*gen-Genchamp-AvgFit"', '3D3+H';...
    '"001-structure/mhn_.*gen-Genchamp-AvgFit"', 'MHN';...
    '"001-structure/mhn-h_.*gen-Genchamp-AvgFit"', 'MHN+H';...
    '"001-structure/mhn3_.*gen-Genchamp-AvgFit"', 'MHN3';... % 15
    '"001-structure/mhn3-h_.*gen-Genchamp-AvgFit"', 'MHN3+H';...
    '"001-structure/rmhn_.*gen-Genchamp-AvgFit"', 'RHN';...
    '"001-structure/rmhn-h_.*gen-Genchamp-AvgFit"', 'RHN+H';...
    '"001-structure/rmhn2_.*gen-Genchamp-AvgFit"', 'RHN2';...
    '"001-structure/rmhn2-h_.*gen-Genchamp-AvgFit"', 'RHN2';... % 20
    '"001-structure/n_.*gen-Genchamp-AvgFit"', 'N';...
});

cd /Users/dk/research/prj/distributed-control/test/var
d2 = quick_multiload({...
    '"001-ini-h.dat"','';...
    '"001-ini.dat"','';...
    '"001-geo1-h.dat"','';...
    '"001-geo1.dat"','';...
    '"001-geo2-h.dat"','';...
    '"001-geo2.dat"','';...
    '"001-3d-h.dat"','';...
    '"001-3d.dat"','';...
    '"001-3d2-h.dat"','';...
    '"001-3d2.dat"','';...
    '"001-3d3-h.dat"','';...
    '"001-3d3.dat"','';...
    '"001-mhn-h.dat"','';...
    '"001-mhn.dat"','';...
    '"001-mhn3-h.dat"','';...
    '"001-mhn3.dat"','';...
    '"001-rmhn-h.dat"','';...
    '"001-rmhn.dat"','';...
    '"001-rmhn2-h.dat"','';...
    '"001-rmhn2.dat"','';...
    '"001-n.dat"','';...
});
cd /Users/dk/research/prj/distributed-control/var

newfigure();
hold on;
for i=1:size(d1,1)
    if(size(d1{i}.data,1) == size(d2{i}.fitness,2))
        x = d2{i}.fitness'./1600;
        y = d1{i}.data(:,end,2)./1600;
        scatter(x,y);
    end 
end
axis([0 1 0 1]);
line([0 1], [0 1]);
xlabel('Sensors disabled');
ylabel('Sensors enabled');
title('Self-organization');
quick_export('../doc/gecco2010/figures/self-org.eps');

%%
% %%
% % Nothing of interest below here...
% 
% %%
% d = quick_plot('data(1,:,1)', 'data(:,:,2)',...
%     {'"009-fewnodes/n-node1_.*gen-Genchamp-AvgFit"', '1';...
%     '"009-fewnodes/n-node2_.*gen-Genchamp-AvgFit"', '2';...
%     '"009-fewnodes/n-node4_.*gen-Genchamp-AvgFit"', '4';...
%     '"009-fewnodes/n-node8_.*gen-Genchamp-AvgFit"', '8';...
%     '"009-fewnodes/n-node16_.*gen-Genchamp-AvgFit"', '16';...
%     '"009-fewnodes/n-node32_.*gen-Genchamp-AvgFit"', '32'},...
%     @(x,y) (y./(2^(x-1)*100)));
% axis([0 200 0 1.05]);
% xlabel('Generations')
% ylabel('Fraction of max. fitness');
% title('NEAT, varying network size');
% quick_export('../doc/gecco2010/figures/009-neat.eps');
% 
% %%
% d = quick_plot('data(1,:,1)', 'data(:,:,2)',...
%     {'"009-fewnodes/hn-node1_.*gen-Genchamp-AvgFit"', '1';...
%     '"009-fewnodes/hn-node2_.*gen-Genchamp-AvgFit"', '2';...
%     '"009-fewnodes/hn-node4_.*gen-Genchamp-AvgFit"', '4';...
%     '"009-fewnodes/hn-node8_.*gen-Genchamp-AvgFit"', '8';...
%     '"009-fewnodes/hn-node16_.*gen-Genchamp-AvgFit"', '16';...
%     '"009-fewnodes/hn-node32_.*gen-Genchamp-AvgFit"', '32'},...
%     @(x,y) (y./(2^(x-1)*100)));
% axis([0 200 0 1.05]);
% xlabel('Generations')
% ylabel('Fraction of max. fitness');
% title('HyperNEAT, varying network size');
% quick_export('../doc/gecco2010/figures/009-hyperneat.eps');
% 
% %%
% d = quick_plot('data(1,:,1)', 'data(:,:,2)',...
%     {'"009-fewnodes/mhn-node1_.*gen-Genchamp-AvgFit"', '1';...
%     '"009-fewnodes/mhn-node2_.*gen-Genchamp-AvgFit"', '2';...
%     '"009-fewnodes/mhn-node4_.*gen-Genchamp-AvgFit"', '4';...
%     '"009-fewnodes/mhn-node8_.*gen-Genchamp-AvgFit"', '8';...
%     '"009-fewnodes/mhn-node16_.*gen-Genchamp-AvgFit"', '16';...
%     '"009-fewnodes/mhn-node32_.*gen-Genchamp-AvgFit"', '32'},...
%     @(x,y) (y./(2^(x-1)*100)));
% axis([0 200 0 1.05]);
% xlabel('Generations')
% ylabel('Fraction of max. fitness');
% title('Multi-agent HyperNEAT, varying network size');
% quick_export('../doc/gecco2010/figures/009-multiagent.eps');
% 
% 
% %%
% d0 = quick_load('"009-fewnodes/hnc-node1_.*gen-Genchamp-AvgFit"');
% d1 = quick_load('"009-fewnodes/hnc-node2_.*gen-Genchamp-AvgFit"');
% d2 = quick_load('"009-fewnodes/hnc-node4_.*gen-Genchamp-AvgFit"');
% d3 = quick_load('"009-fewnodes/hnc-node8_.*gen-Genchamp-AvgFit"');
% d4 = quick_load('"009-fewnodes/hnc-node16_.*gen-Genchamp-AvgFit"');
% d5 = d4;%quick_load('"009-fewnodes/hnc-node32_.*gen-Genchamp-AvgFit"');
% warning('d5==d4/2');
% my_plot6(d0,d1,d2,d3,d4,d5);
% quick_export('../doc/gecco2010/figures/009-hyperneat-c.eps');
% 
% %%
% d0 = quick_load('"009-fewnodes/hni-node1_.*gen-Genchamp-AvgFit"');
% d1 = quick_load('"009-fewnodes/hni-node2_.*gen-Genchamp-AvgFit"');
% d2 = quick_load('"009-fewnodes/hni-node4_.*gen-Genchamp-AvgFit"');
% d3 = quick_load('"009-fewnodes/hni-node8_.*gen-Genchamp-AvgFit"');
% d4 = quick_load('"009-fewnodes/hni-node16_.*gen-Genchamp-AvgFit"');
% d5 = quick_load('"009-fewnodes/hni-node32_.*gen-Genchamp-AvgFit"');
% 
% my_plot6(d0,d1,d2,d3,d4,d5);
% quick_export('../doc/gecco2010/figures/009-hyperneat-id.eps');
% 
% %%
% d0 = quick_load('"010-limited-location/8mod1_.*gen-Genchamp-AvgFit"');
% d1 = quick_load('"010-limited-location/8mod2_.*gen-Genchamp-AvgFit"');
% d2 = quick_load('"010-limited-location/8mod4_.*gen-Genchamp-AvgFit"');
% d3 = quick_load('"010-limited-location/8mod8_.*gen-Genchamp-AvgFit"');
% 
% my_plot4(d0,d1,d2,d3,800);
% quick_export('../doc/gecco2010/figures/010-8nodes.eps');
% 
% %%
% d0 = quick_load('"010-limited-location/8mod1-24_.*gen-Genchamp-AvgFit"');
% d1 = quick_load('"010-limited-location/8mod2-24_.*gen-Genchamp-AvgFit"');
% d2 = quick_load('"010-limited-location/8mod4-24_.*gen-Genchamp-AvgFit"');
% d3 = quick_load('"010-limited-location/8mod8-24_.*gen-Genchamp-AvgFit"');
% 
% my_plot4(d0,d1,d2,d3,800);
% quick_export('../doc/gecco2010/figures/010-8nodes-24.eps');
% 
% %%
% d0 = quick_load('"010-limited-location/16mod1_.*gen-Genchamp-AvgFit"');
% d1 = quick_load('"010-limited-location/16mod2_.*gen-Genchamp-AvgFit"');
% d2 = quick_load('"010-limited-location/16mod4_.*gen-Genchamp-AvgFit"');
% d3 = quick_load('"010-limited-location/16mod8_.*gen-Genchamp-AvgFit"');
% 
% my_plot4(d0,d1,d2,d3,1600);
% quick_export('../doc/gecco2010/figures/010-16nodes.eps');
% 
% %%
% d0 = quick_load('"010-limited-location/32mod1_.*gen-Genchamp-AvgFit"');
% d1 = quick_load('"010-limited-location/32mod2_.*gen-Genchamp-AvgFit"');
% d2 = quick_load('"010-limited-location/32mod4_.*gen-Genchamp-AvgFit"');
% d3 = quick_load('"010-limited-location/32mod8_.*gen-Genchamp-AvgFit"');
% 
% my_plot4(d0,d1,d2,d3,3200);
% quick_export('../doc/gecco2010/figures/010-32nodes.eps');
% 
% 
% %%
% champ_plot('"007-pct-radius/initial_.*gen-Genchamp-AvgFit"', '2d geometry 0');
% %%
% champ_plot('"008-grid/neat_.*gen-Genchamp-AvgFit"', 'neat', 3200);
% %%
% champ_plot('"008-grid/mhn_.*gen-Genchamp-AvgFit"', 'multi-agent', 3200);
% 
% %%
% champ_plot('"008-grid/initial_.*gen-Genchamp-AvgFit"', '2d geometry 0', 3200);
% %%
% champ_plot('"008-grid/ids_.*gen-Genchamp-AvgFit"', 'ids', 3200);
% %%
% 
% d0 = quick_load('"008-grid/geo1_.*gen-Genchamp-AvgFit"'); % 1
% d1 = quick_load('"008-grid/geo2_.*gen-Genchamp-AvgFit"'); % 2
% d2 = quick_load('"008-grid/3d_.*gen-Genchamp-AvgFit"'); % 4
% d3 = quick_load('"008-grid/ids_.*gen-Genchamp-AvgFit"'); % 8
% d4 = quick_load('"008-grid/initial_.*gen-Genchamp-AvgFit"'); % 16
% d5 = quick_load('"008-grid/rxnull_.*gen-Genchamp-AvgFit"'); % 32
% 
% my_plot6b(d0,d1,d2,d3,d4,d5);
% %%
% champ_plot('"008-grid/rxnull_.*gen-Genchamp-AvgFit"', 'rx null', 3200);
% 
% 
% %%
% champ_plot('"009-few-nodes/initial_.*gen-Genchamp-AvgFit"', '2d geometry 0', 800);
% champ_plot('"009-few-nodes/ids_.*gen-Genchamp-AvgFit"', 'ids', 800);
% champ_plot('"009-few-nodes/geo1_.*gen-Genchamp-AvgFit"', '2d geometry 1', 800);
% champ_plot('"009-few-nodes/geo2_.*gen-Genchamp-AvgFit"', '2d geometry 2', 800);
% champ_plot('"009-few-nodes/3d_.*gen-Genchamp-AvgFit"', '3d', 800);
% 
% %%
% champ_plot('"005-substrate-cpl/geo1_.*gen-Genchamp-AvgFit"', '2d geometry 1');
% %%
% champ_plot('"005-substrate-cpl/geo2_.*gen-Genchamp-AvgFit"', '2d geometry 2');
% %%
% champ_plot('"005-substrate-cpl/3d_.*gen-Genchamp-AvgFit"', '3d geometry 1');
% 
% %%
% df = quick_load('"005-substrate-cpl/geo1_.*gen-Genchamp-AvgFit"');
% 
% time = df.data(1,:,1);
% skip = floor(size(time,2)/markers_per_line);
% cf = df.data(:,:,2); % champ fitness
% mcf = mean(cf); % mean champ fitness
% cfse = stderr(cf); % std error of champ fitness
% mpf = mean(df.data(:,:,4)); % mean population fitness
% 
% [fig axes] = newfigure(); box('off'); hold('all');
% 
% p = plot(time, mcf, 'r', time(1:skip:end), mcf(1:skip:end), 'r^',...
%         time, mpf, 'g', time(1:skip:end), mpf(1:skip:end), 'go',...
%         time, mcf+cfse, 'r:', time, mcf-cfse, 'r:',...
%         'MarkerSize', marker_size);
% 
% lh = [p(1) p(3) p(5)];
% legend(lh, 'Champ fit.', 'Mean fit.', 'Std. err.', 'Orientation','horizontal','Location','SouthOutside');
% [lh, oh, ph, s] = legend(axes);
% xlabel('Generation');
% ylabel('Fitness');
% set(oh(5), 'Marker','^');
% set(oh(7), 'Marker','o');
% title('plane');
% 
% 
% %%
% df = quick_load('"001-diameter/min_.*gen-Genchamp-AvgFit"');
% 
% time = df.data(1,:,1);
% skip = floor(size(time,2)/markers_per_line);
% cf = df.data(:,:,2); % champ fitness
% mcf = mean(cf); % mean champ fitness
% cfse = stderr(cf); % std error of champ fitness
% mpf = mean(df.data(:,:,4)); % mean population fitness
% 
% [fig axes] = newfigure(); box('off'); hold('all');
% 
% p = plot(time, mcf, 'r', time(1:skip:end), mcf(1:skip:end), 'r^',...
%         time, mpf, 'g', time(1:skip:end), mpf(1:skip:end), 'go',...
%         time, mcf+cfse, 'r:', time, mcf-cfse, 'r:',...
%         'MarkerSize', marker_size);
% 
% lh = [p(1) p(3) p(5)];
% legend(lh, 'Champ fit.', 'Mean fit.', 'Std. err.', 'Orientation','horizontal','Location','SouthOutside');
% [lh, oh, ph, s] = legend(axes);
% xlabel('Generation');
% ylabel('Fitness');
% set(oh(5), 'Marker','^');
% set(oh(7), 'Marker','o');
% title('plane');
% 
% 
% %%
% d0 = quick_load('"001-diameter/min_.*gen-Genchamp-AvgFit"');
% d1 = quick_load('"001-diameter/max_.*gen-Genchamp-AvgFit"');
% d2 = quick_load('"001-diameter/tgt_.*gen-Genchamp-AvgFit"');
% 
% time = d0.data(1,:,1);
% 
% [fig axes] = newfigure(); box('off'); hold('all');
% 
% p = plot(time, mean(d0.data(:,:,2)),...
%         time, mean(d1.data(:,:,2)),...
%         time, mean(d2.data(:,:,2)));
% 
% legend('min', 'max', 'tgt', 'location', 'eastoutside');
% xlabel('Generation');
% ylabel('Fitness');
% title('diameter');
% 
% %quick_export('../notebook/figures/001-wavesync-fitness.eps');
% 
% %%
% d0 = quick_load('"002-edge/min_.*gen-Genchamp-AvgFit"');
% d1 = quick_load('"002-edge/max_.*gen-Genchamp-AvgFit"');
% d2 = quick_load('"002-edge/tgt_.*gen-Genchamp-AvgFit"');
% 
% time = d0.data(1,:,1);
% 
% [fig axes] = newfigure(); box('off'); hold('all');
% 
% p = plot(time, mean(d0.data(:,:,2)),...
%         time, mean(d1.data(:,:,2)),...
%         time, mean(d2.data(:,:,2)));
% 
% legend('min', 'max', 'tgt', 'location', 'eastoutside');
% xlabel('Generation');
% ylabel('Fitness');
% title('edge');
% 
% %quick_export('../notebook/figures/001-wavesync-fitness.eps');
% 
% %%
% %d0 = quick_load('"001-diameter/min_.*gen-Genchamp-AvgFit"');
% %d1 = quick_load('"001-diameter/max_.*gen-Genchamp-AvgFit"');
% d0 = quick_load('"003-cpl/tgt_.*gen-Genchamp-AvgFit"');
% 
% time = d0.data(1,:,1);
% 
% [fig axes] = newfigure(); box('off'); hold('all');
% 
% p = plot(time, mean(d0.data(:,:,2)));
% 
% legend('tgt', 'location', 'eastoutside');
% xlabel('Generation');
% ylabel('Fitness');
% title('cpl');
% 
% %quick_export('../notebook/figures/001-wavesync-fitness.eps');
% 
% %%
% %d0 = quick_load('"001-diameter/min_.*gen-Genchamp-AvgFit"');
% %d1 = quick_load('"001-diameter/max_.*gen-Genchamp-AvgFit"');
% d0 = quick_load('"004-cc/tgt_.*gen-Genchamp-AvgFit"');
% 
% time = d0.data(1,:,1);
% 
% [fig axes] = newfigure(); box('off'); hold('all');
% 
% p = plot(time, mean(d0.data(:,:,2)));
% 
% legend('tgt', 'location', 'eastoutside');
% xlabel('Generation');
% ylabel('Fitness');
% title('cc');
% 
% %quick_export('../notebook/figures/001-wavesync-fitness.eps');
% 
% 
% %%
% df = quick_load('"001-plane/manet.*gen-Genchamp-AvgFit"');
% 
% time = df.data(1,:,1);
% skip = floor(size(time,2)/markers_per_line);
% cf = df.data(:,:,2); % champ fitness
% mcf = mean(cf); % mean champ fitness
% cfse = stderr(cf); % std error of champ fitness
% mpf = mean(df.data(:,:,4)); % mean population fitness
% 
% [fig axes] = newfigure(); box('off'); hold('all');
% 
% p = plot(time, mcf, 'r', time(1:skip:end), mcf(1:skip:end), 'r^',...
%         time, mpf, 'g', time(1:skip:end), mpf(1:skip:end), 'go',...
%         time, mcf+cfse, 'r:', time, mcf-cfse, 'r:',...
%         'MarkerSize', marker_size);
% 
% lh = [p(1) p(3) p(5)];
% legend(lh, 'Champ fit.', 'Mean fit.', 'Std. err.', 'Orientation','horizontal','Location','SouthOutside');
% [lh, oh, ph, s] = legend(axes);
% xlabel('Generation');
% ylabel('Fitness');
% set(oh(5), 'Marker','^');
% set(oh(7), 'Marker','o');
% title('plane');
% 
% %quick_export('../notebook/figures/001-wavesync-fitness.eps');
% %%
% %%
% dir='002-small';
% d0 = quick_load([ '"' dir '/manet16_.*gen-Genchamp-AvgFit"']);
% d1 = quick_load([ '"' dir '/manet8_.*gen-Genchamp-AvgFit"']);
% d2 = quick_load([ '"' dir '/manet4_.*gen-Genchamp-AvgFit"']);
% d3 = quick_load([ '"' dir '/manet3_.*gen-Genchamp-AvgFit"']);
% 
% [fig axes] = newfigure(); box('off'); hold('all');
% time = d0.data(1,:,1);
% 
% p = plot(time, (mean(d0.data(:,:,4))),...
%         time, (mean(d1.data(:,:,4))),...
%         time, (mean(d2.data(:,:,4))),...
%         time, (mean(d3.data(:,:,4))));
%     
% legend('0.4', '0.2', '0.1', '0.05',...
%     'location', 'eastoutside');
% ylabel('Fitness');
% xlabel('Time (updates)');
% title('small networks');
% 
% 
% 
% %%
% dir='003-kill';
% d0 = quick_load([ '"' dir '/manet16_.*gen-Genchamp-AvgFit"']);
% d1 = quick_load([ '"' dir '/manet8_.*gen-Genchamp-AvgFit"']);
% d2 = quick_load([ '"' dir '/manet4_.*gen-Genchamp-AvgFit"']);
% d3 = quick_load([ '"' dir '/manet3_.*gen-Genchamp-AvgFit"']);
% 
% [fig axes] = newfigure(); box('off'); hold('all');
% time = d0.data(1,:,1);
% 
% p = plot(time, (mean(d0.data(:,:,4))),...
%         time, (mean(d1.data(:,:,4))),...
%         time, (mean(d2.data(:,:,4))),...
%         time, (mean(d3.data(:,:,4))));
%     
% legend('0.4', '0.2', '0.1', '0.05',...
%     'location', 'eastoutside');
% ylabel('Fitness');
% xlabel('Time (updates)');
% title('kill');
% 
% 
% 
% %%
% dir='004-degree';
% d0 = quick_load([ '"' dir '/manet_.*gen-Genchamp-AvgFit"']);
% 
% [fig axes] = newfigure(); box('off'); hold('all');
% time = d0.data(1,:,1);
% 
% p = plot(time, (mean(d0.data(:,:,4))));
%     
% legend('degree 3',...
%     'location', 'eastoutside');
% ylabel('Fitness');
% xlabel('Time (updates)');
% title('');
% 
% %%
% dir='005-bb';
% d0 = quick_load([ '"' dir '/manet_.*gen-Genchamp-AvgFit"']);
% 
% [fig axes] = newfigure(); box('off'); hold('all');
% time = d0.data(1,:,1);
% 
% p = plot(time, (mean(d0.data(:,:,4))));
%     
% legend('degree 3',...
%     'location', 'eastoutside');
% ylabel('Fitness');
% xlabel('Time (updates)');
% title('');