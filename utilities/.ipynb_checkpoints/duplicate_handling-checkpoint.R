# library
library(annotables)
library(dplyr)

# map the id
get_duplicated_idx <- function(nums) {
        nums = sort(nums)
        slow = 1
        duplicated_idx = c()
        for (fast in 2:length(nums)){
            if (nums[fast] != nums[slow]) {
                slow = slow + 1
                nums[slow] = nums[fast]
                }  else {
            duplicated_idx <- c(duplicated_idx, fast)
            }
        }
        #return(slow + 1)
        return(duplicated_idx)
}
                          
        
    
    
