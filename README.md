# The Moonbounce Project

<!--- These are examples. See https://shields.io for others or to customize this set of shields. You might want to include dependencies, project status and licence info here --->

![GitHub repo size](https://img.shields.io/github/repo-size/Alek99/Moonbounce)
![GitHub contributors](https://img.shields.io/github/contributors/Alek99/Moonbounce)
![GitHub stars](https://img.shields.io/github/stars/Alek99/Moonbounce?style=social)
![GitHub forks](https://img.shields.io/github/forks/Alek99/Moonbounce?style=social)

## Project Overview
We propose to use the Green Bank Telescope (GBT) located in West Virginia, to observe Earth’s radio leakage radiation as reflected from the moon in the 290-1230 MHz range quarterly and yearly. These observations will give an indication of how the characteristics and detectability of Earth’s radio leakage, as seen by an external observer, have changed since a similar measurement was last performed in 2012 and before that in 1978. The Moonbounce project has previously been awarded observation time on the Green Bank Telescope in the fall of 2018 but observations did not come to fruition due to the priority level of the awarded time. The team has since re-submitted the proposal in January of 2019. 

## Results
Given a specific observation time a corresponding area on earth is able to emmit signals that can bounce off the moon and be collected at our reciever (GBT) which is locate in West Virginia. The created program takes in an obervation time and outputs an area that the possible collected signal could have originated from. In the generated graph below the lighter region represents the designated area of possible signals.

![alt text](https://github.com/Alek99/Moonbounce/blob/master/October6ObservationExample.png)

## WebApp

The WebApp displays the same functionality listed above but in a dynamic way. Once running simply input an observation time and a graph of moon-earth visibility will be displayed. Additionally another feature currently in progress is a database of our observation times which can be filtered and sorted based on criteria such as observation date, length, and desired signal range. Example of webapp display below where the yellow represents the moon's overhead path and red and yellow lines represent where visibilty starts and ends.

![alt text](https://github.com/Alek99/Moonbounce/blob/master/Screen%20Shot%202020-04-15%20at%2011.54.22%20AM.png)


## Prerequisites

Before you begin, ensure you have met the following requirements:
<!--- These are just example requirements. Add, duplicate or remove as required --->
* You have installed the latest version of Flask
* You have a Windows/Linux/Mac machine.

## Installing and Running Moonbounce WebApp

To install the Moonbounce Flask app, once in the Moonbounce WebApp folder follow these steps:

$ pip3 install -r requirements.txt 

Set the FLASK_APP environment variable 

* $ export FLASK_APP=run.py 

Start the application (development mode) 

* $ flask run 

Access the dashboard in browser: http://127.0.0.1:5000/ 
Create an account and log in 

## Contributors

Thanks to the following people who have contributed to this project:

* [@scottydocs](https://github.com/scottydocs) 📖
* [@cainwatson](https://github.com/cainwatson) 🐛
* [@calchuchesta](https://github.com/calchuchesta) 🐛

You might want to consider using something like the [All Contributors](https://github.com/all-contributors/all-contributors) specification and its [emoji key](https://allcontributors.org/docs/en/emoji-key).

## Contact

If you want to contact me you can reach me at 17petuskey@berkeley.edu.

## License
<!--- If you're not sure which open license to use see https://choosealicense.com/--->

This project uses the following license: [<license_name>](<link>).
