function T = tensor3_to_tensor22(H)

T(:,:,1,1) = H(:,:,1);
T(:,:,2,2) = H(:,:,2);
T(:,:,1,2) = H(:,:,3);
T(:,:,2,1) = H(:,:,3);
    
end
