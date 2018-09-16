# Performance Lasertools

## Settings
- File: drawing_test_Panther_small.svg
- 4 core i5 Laptop

## Test 1
### Settings
- Beam with: 0.5
- Passes: 1
- Accuracy: 20
- Create log: true
- Multithreading: false

### Results
- 16.3408961296 s for infill 
- 40.6142289639 s for parameters 
    - 40.1177051067 s for finishing path 

## Test 2
### Settings
- Beam with: 0.5
- Passes: 1
- Accuracy: 20
- Create log: true
- Multithreading: true

### Results
- 5.30531001091 s for infill 
- 39.2469501495 s for parameters -- not multithreaded yet!
    - 38.9640929699 s for finishing path

## Test 3
### Settings
- Beam with: 0.5
- Passes: 1
- Accuracy: 10
- Create log: true
- Multithreading: false

### Results
- 16.3090660572 s for infill 
- 34.0203380585 s for parameters 

## Test 4
### Settings
- Beam with: 0.5
- Passes: 1
- Accuracy: 2
- Create log: true
- Multithreading: false

### Results
- 16.368336916 s for infill 
- 33.019288063 s for parameters 

## Test 5
### Settings
- Beam with: 0.5
- Passes: 1
- Accuracy: 2
- Create log: true
- Multithreading: true
- biarc-max-split-depth: 1

### Results
- 5.28588199615 s for infill 
- 38.9247260094 s for parameters 

## Test 6
### Settings
- Beam with: 0.5
- Passes: 1
- Accuracy: 2
- Create log: true
- Multithreading: true
- biarc-max-split-depth: 100

### Results
- 5.22256612778 s for infill 
- 32.8186120987 s for parameters 

## Test 7
### Settings
- Beam with: 0.5
- Passes: 1
- Accuracy: 2
- Create log: true
- Multithreading: true
- engraving accuracy: 0.2

### Results
- 5.63703203201 s for infill 
- 15.0216479301 s for parameters 

# Test Removing emty G00 Instructions


## Test With empty G00 Instructions
### Settings
- Beam with: 0.35
- Passes: 1
- Create log: true
- Multithreading: true

### Results
- 11695 Lines of G-Code
- 8.39188909531 s for infill
- 19.0251300335 s for parameters

## Test Without empty G00 Instructions
### Settings
- Beam with: 0.35
- Passes: 1
- Create log: true
- Multithreading: true

### Results
- 10247 Lines of G-Code
- 8.09686088562 s for infill
- 18.9610948563 s for parameters 

## Conclusion
no significant perfomance impact; generates cleaner code

