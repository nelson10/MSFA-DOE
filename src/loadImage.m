function [I] = loadImage(detector,l,dynamicRange)
if(detector <= 2)
    alldataset = {'balloons_ms','beads_ms','cd_ms','chart','clay_ms','cloth_ms','egyptian_statue_ms','feathers_ms','flowers_ms','glass_tiles_ms','pompoms_ms','sponges_ms','stuffed_toys_ms','superballs_ms','thread_spools_ms','fake_and_real_beers_ms','face_ms','real_and_fake_peppers_ms','real_and_fake_apples_ms','photo_and_face_ms','paints_ms','oil_painting_ms','jelly_beans_ms','hairs_ms','fake_and_real_tomatoes_ms','fake_and_real_sushi_ms','fake_and_real_strawberries_ms','fake_and_real_peppers_ms','fake_and_real_lemons_ms','fake_and_real_lemon_slices_ms','fake_and_real_food_ms','watercolors_ms'};
    dataset = alldataset{l};
    load(dataset);
    I = mat2gray(double(hyperimg))*dynamicRange;
elseif(detector == 3)
    if(l<10)
        I = imread("kodim0"+num2str(l)+".png"); % get linear RGB image
    else
        I = imread("kodim"+num2str(l)+".png"); % get linear RGB image
    end
    I = mat2gray(double(I))*dynamicRange;
end
end