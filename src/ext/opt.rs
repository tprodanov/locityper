

/// Compile time zero-sized structure, used as a more readable alternative to const bool generic.
pub trait Flag : Copy + Default {
    fn enabled(self) -> bool;
}

#[derive(Clone, Copy, Default)]
pub struct Off;

impl Flag for Off {
    #[inline(always)]
    fn enabled(self) -> bool { false }
}

#[derive(Clone, Copy, Default)]
pub struct On;

impl Flag for On {
    #[inline(always)]
    fn enabled(self) -> bool { true }
}

/// Compile-time option.
pub trait Optional<T> : Copy {
    type Flag : Flag;

    #[inline(always)]
    fn flag(self) -> Self::Flag {
        Self::Flag::default()
    }

    #[inline(always)]
    fn enabled(self) -> bool {
        self.flag().enabled()
    }

    fn get(self) -> Option<T>;
}

impl<T> Optional<T> for Off {
    type Flag = Off;

    #[inline(always)]
    fn get(self) -> Option<T> { None }
}

#[derive(Clone, Copy)]
pub struct With<T: Copy>(pub T);

impl<T: Copy> Optional<T> for With<T> {
    type Flag = On;

    #[inline(always)]
    fn get(self) -> Option<T> { Some(self.0) }
}
