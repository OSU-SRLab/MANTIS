#pragma once

#include <string>

using std::string;

namespace mantis
{
	class CircularLinkedList
	{
		private:
		class Node
		{
			private:
			char element;
			Node* next;
			
			public:
			void setElement(const char& c);
			char getElement() const;
			void setNext(Node* n);
			Node* getNext() const;
		};

		int k;
		unsigned int curSeqLength;
		int numIgnoreChars;
		string ignoreChar;
		unsigned int lastSeqLength;
		int repeatCount;
		Node* first;
		Node* cur;
		string kmer;
		int minBases;
		int maxBases;
		int minRepeats;
		
		bool isCharInIgnore(const char& c) const;
		bool isCharN(const char& c) const; //much faster than isCharInIgnore if we're only ignoring N's
		bool seqHasIgnoreChrs() const;
		void pullKmer();
		
		public:

		/*
		* Makes a new Circular Linked List for finding repeat patterns.
		* @warning user_ignoreChar is currently ignored, and any repeat patterns containing N will not be reported.
		* @param user_k the value of k (repeat sequence size) for this CLL to detect
		* @param user_minBases the minimum number of bp a repeat sequence must cover to be reported
		* @param user_maxBases the maximum number of bp a repeat sequence must cover to be reported
		* @param user_minRepeats the minimum number of times a repeat sequence must repeat to be reported
		* @param user_ignoreChar a string of characters, the presence of which in a sequence will cause it not
		*        to be reported. Currently IGNORED.
		*/
		CircularLinkedList(int user_k, int user_minBases, int user_maxBases, int user_minRepeats, const string& user_ignoreChar);

		~CircularLinkedList();

		/*
		* Considers the next element in the sequence being scanned, and reports if a sequence matching criteria has
		*        just been ended. If true, then use getLastLength() and getKmer() to report it.
		* @param c the next char in the sequence
		* @return true if a sequence matching criteria was just ended, false otherwise
		*/
		bool next(const char& c);

		/*
		* Updates this CLL with the next element in the sequence being scanned, but does NOT check if a sequence
		*        matching criteria has just been ended.
		* @param c the next char in the sequence
		*/
		void nextNoResultCheck(const char& c);

		/*
		* Reports the value of k this CLL was created with.
		* @return k
		*/
		int getK() const;

		/*
		* Causes the CLL to check if the current sequence being considered meets criteria. Call when at the end
		*        of the sequence being scanned, to ensure that repeat sequences continuing to the end are properly
		*        detected.
		* @return true if a sequence matching criteria was just ended, false otherwise
		*/
		bool endSequence(); //exposed so that the user can let it know when it's out of data

		/*
		* Reports the length in bp of the last repeat sequence matching criteria found. Only meaningful if next()
		*        or endSequence() returned true.
		* @return length of last repeat sequence
		*/
		unsigned int getLastLength() const;

		/*
		* Returns the string of characters that this CLL is ignoring sequences with.
		* @return ignore characters string
		*/
		string getIgnoreChars() const;

		/*
		* Changes the string of characters that this CLL is ignoring sequences with.
		* @param user_ignoreChar the new ignore characters string
		*/
		void setIgnoreChars(const string& user_ignoreChar);

		/*
		* Reports the repeat sequence of the last repeat sequence matching criteria found. Only meaningful if
		*        next() or endSequence() returned true.
		* @return the repeat sequence
		*/
		string getKmer() const;

		/*
		* Reports the number of times the last repeat sequence matching criteria repeated. Only meaningful if
		*        next() or endSequence() returned true.
		* @return the repeat count
		*/
		int getRepeatCount() const;
	}; //end class CircularLinkedList

	inline void CircularLinkedList::Node::setElement(const char& c)
	{
		element = c;
	}

	inline char CircularLinkedList::Node::getElement() const
	{
		return element;
	}

	inline void CircularLinkedList::Node::setNext(CircularLinkedList::Node* n)
	{
		next = n;
	}

	inline CircularLinkedList::Node* CircularLinkedList::Node::getNext() const
	{
		return next;
	}

	inline bool CircularLinkedList::isCharN(const char& c) const
	{
		return (c == 'N');
	}
} //end namespace mantis
