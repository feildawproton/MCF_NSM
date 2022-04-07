# MCF_NSM
Minimum cost flow using network simplex method. 2018 for MSIM722.

This is an old school project that I decided to upload for safe keeping.

## Minimum Cost Flow without Upper and Lower Bounds
- Implemented the network simplex method for unbounded minimum cost flow problems
- Tested on 3 small and 1 medium size networks
- Many of the larger problems were in a custom format and I did not have the time to test them
- Output either validated from source or by constructed solution
- Troubled development would indicate that this program is not robust
  - Many manual edge case checks
  - Not nearly as well planned, organized, or executed as my GA program  
- Tedious input method

## Many checks and helper functionsd during the calculation of w

## Creating a new basis
- Doing so in an automated manner requires checking
  - Checking which nodes are dead end nodes
  - Checking which links are connected to dead end nodes (static links)
    - Ensuring that static link flows will be satisfied
      - Based on their connected dead end node and the other link flows
  - Updating link flows from know previous basis flows

## example 1:
![image](https://user-images.githubusercontent.com/56926839/162252777-723f4e4e-83e0-4ee7-b91e-754351a8df0a.png)

## example 2:
![image](https://user-images.githubusercontent.com/56926839/162252816-68859273-4ae2-4444-abee-ab31869f89f2.png)

## example 3 from sas.com:
![image](https://user-images.githubusercontent.com/56926839/162252864-32d29e54-abf6-4dd4-94e5-9ef793986522.png)

## example 4: Scaling test with 33 nodes, 51 links and a known solutions
![image](https://user-images.githubusercontent.com/56926839/162252969-709ea0e2-a5bf-44f7-a6e6-f56a1ff6e842.png)




