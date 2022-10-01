#include <memory>
#include <iostream>
#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/thread/thread.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <SFML/Graphics.hpp>

#include "../environment/ping-pong-environment.h"

#define RASTER_HALFSIZE 100   // It is for the range [0, 1] while the whole picture range is [-1, 1].
#define SPOT_HALFSIZE_PIXEL 10

using namespace std;
using namespace boost::interprocess;

class EnvironmentState
{
	unique_ptr<shared_memory_object> shm;
	unique_ptr<mapped_region>        region;
public:
	volatile pair<float, float> *pprr_Ball;
	volatile float              *prRacket = NULL;
	EnvironmentState(bool bDestroyOldState)
	{
		if (!bDestroyOldState)
			do {
				try {
					//Create a shared memory object.
					shm.reset(new shared_memory_object(open_only, ENVIRONMENT_STATE_SHARED_MEMORY_NAME, read_only));

					//Map the whole shared memory in this process
					region.reset(new mapped_region(*shm, read_only));
					pprr_Ball = (volatile pair<float, float> *)region->get_address();
					prRacket = &((volatile pair<pair<float, float>, float> *)pprr_Ball)->second;
				} catch (...) {
					boost::this_thread::sleep(boost::posix_time::seconds(1));
				}
			} while (!prRacket);
		else {
			string strSharedMemoryName = ENVIRONMENT_STATE_SHARED_MEMORY_NAME;
			bool bExists = true;
			do {
				try {
					//Create a shared memory object.
					shm.reset(new shared_memory_object(open_only, strSharedMemoryName.c_str(), read_only));
					shared_memory_object::remove(strSharedMemoryName.c_str());
					++strSharedMemoryName.front();
				} catch (...) {
					bExists = false;
				}
			} while (bExists);
		}
	}
	~EnvironmentState() {}
};

inline float rPixelX(float rPhysical) {return (rPhysical + 1.F) * RASTER_HALFSIZE;}
inline float rPixelY(float rPhysical) {return (rPhysical + 1.F) * RASTER_HALFSIZE;}

int main(int ARGC, char *ARGV[])
{
	bool bDestroyOldState = ARGC == 2;
	cout << "Waiting for ping-pong data availability...\n";
	EnvironmentState es(bDestroyOldState);
	if (bDestroyOldState) {
		cout << "All resident ping-pong data are destroyed\n";
		exit(0);
	}
	cout << "ping-pong data are obtained\n";
	// Create the main window
	sf::RenderWindow window(sf::VideoMode(RASTER_HALFSIZE * 2, RASTER_HALFSIZE * 2), "ping-pong", sf::Style::Close);
	sf::CircleShape circle;
	circle.setRadius(SPOT_HALFSIZE_PIXEL);
	circle.setFillColor(sf::Color::White);
	circle.setOrigin(SPOT_HALFSIZE_PIXEL, SPOT_HALFSIZE_PIXEL);
	sf::RectangleShape rectangle;
	rectangle.setSize(sf::Vector2f(3, RACKET_SIZE * RASTER_HALFSIZE));
	rectangle.setOrigin(0, RACKET_SIZE / 2 * RASTER_HALFSIZE);
	rectangle.setFillColor(sf::Color::White);

	// Start the game loop
	while (window.isOpen())
	{
		// Process events
		sf::Event event;
		while (window.pollEvent(event))
		{
			// Close window: exit
			if (event.type == sf::Event::Closed)
				window.close();
		}

		// Clear screen
		window.clear();
		circle.setPosition(rPixelX(es.pprr_Ball->first), rPixelY(es.pprr_Ball->second));
		window.draw(circle);
		rectangle.setPosition(0, rPixelY(*es.prRacket));
		window.draw(rectangle);

		// Update the window
		window.display();
	}
	return EXIT_SUCCESS;
}
