
function weights = stepWeights(total, step)
    weights = ones(total, 1);
    for i = 1:total
        if i < step
            weights(i) = weights(i) - ((1-0.2)/step)*i;
        else
            weights(i) = 0.2;
        end
    end
end