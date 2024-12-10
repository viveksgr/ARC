function projV1onV2 = ARC_projectV1onV2(V1, V2)
    % Ensure V1 and V2 are the same size
    if length(V1) ~= length(V2)
        error('V1 and V2 must be of the same size');
    end

    % Calculate the projection of V1 onto V2
    projV1onV2 = (dot(V1, V2) / dot(V2, V2)) * V2;
end
