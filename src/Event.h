#ifndef _EVENT_H_
#define _EVENT_H_

#include <vector>
#include <tuple>
#include <cstddef>

// Event type enum
enum event {
  INSERTION,
  DELETION
};

// Event structure representing a single indel event
struct Event {
    event type;        // INSERTION or DELETION
    size_t position;   // Position in current sequence
    size_t length;     // Length of insertion/deletion
};

// Sequence of events along a branch
typedef std::vector<Event> EventSequence;

// Map from node ID to (event sequence, final sequence size)
// Using vector since nodes are numbered 0 to n-1
typedef std::vector<EventSequence> EventMap;

#endif // _EVENT_H_
