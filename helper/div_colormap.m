%% SPECIFY DIVERGING COLORMAP
%indices = linspace(0,1,26);
%rgb1 = [0.230, 0.299, 0.754];
%rgb2 = [0.706, 0.016, 0.150];
%div_colormap = diverging_map(indices,rgb1,rgb2);

function map = div_colormap 
map = [ ...
    0.2300    0.2990    0.7540; ...
    0.2769    0.3684    0.8144; ...
    0.3258    0.4362    0.8675; ...
    0.3768    0.5016    0.9125; ...
    0.4297    0.5640    0.9488; ...
    0.4841    0.6225    0.9756; ...
    0.5393    0.6764    0.9927; ...
    0.5948    0.7248    0.9998; ...
    0.6496    0.7671    0.9969; ...
    0.7029    0.8026    0.9839; ...
    0.7539    0.8306    0.9613; ...
    0.8016    0.8508    0.9294; ...
    0.8453    0.8627    0.8888; ...
    0.8869    0.8555    0.8372; ...
    0.9223    0.8293    0.7784; ...
    0.9480    0.7951    0.7172; ...
    0.9638    0.7533    0.6546; ...
    0.9701    0.7043    0.5914; ...
    0.9668    0.6487    0.5285; ...
    0.9543    0.5870    0.4666; ...
    0.9328    0.5197    0.4064; ...
    0.9026    0.4469    0.3486; ...
    0.8643    0.3684    0.2936; ...
    0.8183    0.2826    0.2419; ...
    0.7653    0.1831    0.1939; ...
    0.7060    0.0161    0.1500];
end